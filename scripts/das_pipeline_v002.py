import argparse
import glob
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time

# m_input_dir = "/people/whit040/Passlid/input"
# m_output_dir = "/people/whit040/Passlid/output"
# m_config_file = "/people/whit040/Passlid/das_config.txt"
# m_param_file = "/people/whit040/Passlid/das_param.txt"

# m_input_dir = None
# m_output_dir = None
# m_config_file = None
# m_param_file = None

m_config = {}
m_param = {}
m_read_pairs = {}


def init_run_info():
    init_read_pairs()
    init_config()
    init_params()


def init_read_pairs():
    if not os.path.exists(m_input_dir):
        raise Exception("The input directory does not exist.")

    if not os.path.exists(m_output_dir):
        os.makedirs(m_output_dir)

    forward_files = glob.glob(m_input_dir + "/*_R1.fastq")
    for forward_file in forward_files:
        base_name = os.path.basename(forward_file)
        base_name = base_name.replace("_R1.fastq", "")

        reverse_file = m_input_dir + "/" + base_name + "_R2.fastq"
        if not os.path.exists(reverse_file):
            raise Exception("The reverse file does not exist.")

        m_read_pairs[base_name] = (forward_file, reverse_file)


def init_config():
    if not os.path.exists(m_config_file):
        raise Exception("Configuration file does not exist.")

    fIn = open(m_config_file)

    for line in fIn:
        line = line.strip()
        if re.match("#", line) or line == "":
            continue

        line_data = line.split(":")
        key = line_data[0]
        key = key.strip()

        val = line_data[1]
        val = val.strip()

        m_config[key] = val

    fIn.close()

    m_config['INPUT_DIR'] = m_input_dir
    m_config['OUTPUT_DIR'] = m_output_dir


def init_params():
    if not os.path.exists(m_param_file):
        raise Exception("Parameter file does not exist.")

    fIn = open(m_param_file)

    for line in fIn:
        line = line.strip()
        if re.match("#", line) or line == "":
            continue

        line_data = line.split(":")
        program_name = line_data[0]
        if not program_name in m_param.keys():
            m_param[program_name] = {}

        line_data = line_data[1].split(" ")
        param_name = line_data[0]
        param_val = line_data[1:]
        param_val = " ".join(param_val)

        m_param[program_name][param_name] = param_val

    fIn.close()


def merge_reads(read_pair_id):
    flash_dir = get_flash_dir(read_pair_id)

    r1_file = get_forward_read_path(read_pair_id)
    r2_file = get_reverse_read_path(read_pair_id)

    output_prefix = read_pair_id

    if os.path.exists(flash_dir):
        shutil.rmtree(flash_dir)

    os.makedirs(flash_dir)

    start_time = time.time()

    the_args = [m_config['FLASH_EXECUTABLE'], r1_file, r2_file, "-m", m_param['flash']['min_overlap'], "-M",
                m_param['flash']['max_overlap'], "-x", m_param['flash']['mismatch_ratio'],
                "-p", m_param['flash']['phred_offset'], "-o", output_prefix, "-d", flash_dir]

    the_cmd = " ".join(the_args)

    log_file_name = read_pair_id + "_Flash.log"
    log_file_path = os.path.join(flash_dir, log_file_name)
    flog = open(log_file_path, 'w')
    subprocess.call(the_cmd, shell=True, stdout=flog)
    flog.close()

    end_time = time.time()
    # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


def contaminant_db_exists(db_name):
    db_file_post_fixes = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    db_exists = True
    for post_fix in db_file_post_fixes:
        file_name = db_name + post_fix
        file_path = os.path.join(m_config['CONTAMINANT_DB_DIR'], file_name)

        if not os.path.exists(file_path):
            db_exists = False
            continue

    return db_exists


def build_contaminant_db(db_name):
    start_time = time.time()

    fasta_file_name = db_name + ".fa"
    fasta_file_path = os.path.join(m_config['CONTAMINANT_DB_DIR'], fasta_file_name)

    if not os.path.exists(fasta_file_path):
        raise Exception("FASTA file for the contaminant database does not exist.")

    bt2_index_base = os.path.join(m_config['CONTAMINANT_DB_DIR'], db_name)
    
    the_args = [m_config['BW2_BUILD_EXECUTABLE'], fasta_file_path, bt2_index_base]
    # the_cmd = " ".join(the_args)
    subprocess.call(the_args, shell=False)

    end_time = time.time()
    # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


def find_contaminants(read_pair_id, db_name, step_id, final_step):
    # print("Searching contaminant database using Bowtie2 with " + read_pair_id + ", " + db_name + " ..")

    start_time = time.time()

    num_threads = m_param['flash']['num_threads']

    contaminant_db_index_base = os.path.join(m_config['CONTAMINANT_DB_DIR'], db_name)
    flash_dir = get_flash_dir(read_pair_id)
    decon_dir = get_decon_dir(read_pair_id)

    #-------------------------------------------------------------------------------------------------------------------
    # Create or open the  log file
    #-------------------------------------------------------------------------------------------------------------------

    log_file_name = read_pair_id + "_Decon.log"
    log_file_path = os.path.join(decon_dir, log_file_name)

    if not os.path.exists(log_file_path):
        separator = "=" * 100

        flog = open(log_file_path, 'w')
        flog.write(separator + "\n")
        flog.write("Read pair = " + read_pair_id + "\n")
        flog.write(separator + "\n\n")

        flog.flush()
    else:
        flog = open(log_file_path, 'a')

    #-------------------------------------------------------------------------------------------------------------------
    # Find contaminants in extended fragments
    #-------------------------------------------------------------------------------------------------------------------

    separator = "-" * 100

    decon_info = []
    decon_info.append(["Extended_Fragments", "extendedFrags", "Extended_Frags"])
    decon_info.append(["Not Combined 1", "notCombined_1", "Not_Combined_1"])
    decon_info.append(["Not Combined 2", "notCombined_2", "Not_Combined_2"])

    for base_names in decon_info:
        if step_id == 1:
            bw2_input_file_name = read_pair_id + "." + base_names[1] + ".fastq"
            bw2_input_file_path = os.path.join(flash_dir, bw2_input_file_name)
        else:
            prev_step_id = step_id - 1
            prev_step_id = str(prev_step_id)

            the_file_pattern = "Step_" + prev_step_id.zfill(3) + "*_Unaligned_" + base_names[2] + ".fastq"
            the_file_pattern = os.path.join(decon_dir, the_file_pattern)
            bw2_input_file_path = glob.glob(the_file_pattern)
            bw2_input_file_path = bw2_input_file_path[0]

        curr_step_id = str(step_id)

        if final_step:
            step_prefix = "Final_"
        else:
            step_prefix = "Step_" + curr_step_id.zfill(3) + "_"

        bw2_output_file_name = step_prefix + read_pair_id + "_" + db_name + "_" + base_names[2] + ".sam"
        bw2_output_file_path = os.path.join(decon_dir, bw2_output_file_name)

        flog.write(separator + "\n")
        flog.write("Search :: Contaminant DB = " + db_name + ", " + base_names[0] + "\n")
        flog.write(separator + "\n\n")

        flog.flush()

        the_args = [m_config['BW2_EXECUTABLE'], "-p", num_threads, "-x", contaminant_db_index_base, "-q", bw2_input_file_path, "-S", bw2_output_file_path, "--very-sensitive-local"]
        the_cmd = " ".join(the_args)

        subprocess.call(the_cmd, shell=True, stderr=flog)

        flog.write("\n")
        flog.flush()

        flog.write(separator + "\n")
        flog.write("Removal :: Contaminant DB = " + db_name + ", " + base_names[0] + "\n")
        flog.write(separator + "\n\n")

        flog.flush()

        picard_sam_file_name  =step_prefix + read_pair_id + "_" + db_name + "_Unaligned_" + base_names[2] + ".sam"
        picard_sam_file_path = os.path.join(decon_dir, picard_sam_file_name)

        the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Unaligned", "I=" + bw2_output_file_path, ">", picard_sam_file_path]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True,  stderr=flog)

        picard_fastq_file_name = step_prefix + read_pair_id + "_" + db_name + "_Unaligned_" + base_names[2] + ".fastq"
        picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

        the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True,  stderr=flog)

        flog.write("\n")

    # #-------------------------------------------------------------------------------------------------------------------
    # # Find contaminants in non-extended forward (R1) fragments
    # #-------------------------------------------------------------------------------------------------------------------
    #
    # flash_file_name = read_pair_id + ".notCombined_1.fastq"
    # flash_file_path = os.path.join(flash_dir, flash_file_name)
    #
    # sam_file_name = read_pair_id + "_" + db_name + "_Not_Combined_1.sam"
    # sam_file_path = os.path.join(decon_dir, sam_file_name)
    #
    # flog.write(separator + "\n")
    # flog.write("Contaminant DB = " + db_name + ", Not Combined 1\n")
    # flog.write(separator + "\n\n")
    #
    # flog.flush()
    #
    # the_args = [m_config['BW2_EXECUTABLE'], "-p", num_threads, "-x", contaminant_db_index_base, "-q", flash_file_path, "-S", sam_file_path, "--very-sensitive-local"]
    # the_cmd = " ".join(the_args)
    #
    # subprocess.call(the_cmd, shell=True, stderr=flog)
    # flog.write("\n")
    #
    # #-------------------------------------------------------------------------------------------------------------------
    # # Find contaminants from non-extended reverse (R2) fragments
    # #-------------------------------------------------------------------------------------------------------------------
    #
    # flash_file_name = read_pair_id + ".notCombined_2.fastq"
    # flash_file_path = os.path.join(flash_dir, flash_file_name)
    #
    # sam_file_name = read_pair_id + "_" + db_name + "_Not_Combined_2.sam"
    # sam_file_path = os.path.join(decon_dir, sam_file_name)
    #
    # flog.write(separator + "\n")
    # flog.write("Contaminant DB = " + db_name + ", Not Combined 2\n")
    # flog.write(separator + "\n\n")
    #
    # flog.flush()
    #
    # the_args = [m_config['BW2_EXECUTABLE'], "-p", num_threads, "-x", contaminant_db_index_base, "-q", flash_file_path, "-S", sam_file_path, "--very-sensitive-local"]
    # the_cmd = " ".join(the_args)
    #
    # subprocess.call(the_cmd, shell=True, stderr=flog)
    # flog.write("\n")
    #
    # #-------------------------------------------------------------------------------------------------------------------

    flog.close()

    end_time = time.time()
    # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


# def find_contaminants(read_pair_id, db_name):
#     # print("Searching contaminant database using Bowtie2 with " + read_pair_id + ", " + db_name + " ..")
#
#     start_time = time.time()
#
#     num_threads = m_param['flash']['num_threads']
#
#     contaminant_db_index_base = os.path.join(m_config['CONTAMINANT_DB_DIR'], db_name)
#     flash_dir = get_flash_dir(read_pair_id)
#     decon_dir = get_decon_dir(read_pair_id)
#
#     #-------------------------------------------------------------------------------------------------------------------
#     # Create or open the  log file
#     #-------------------------------------------------------------------------------------------------------------------
#
#     log_file_name = read_pair_id + "_Decon_Contaminant_Search.log"
#     log_file_path = os.path.join(decon_dir, log_file_name)
#
#     if not os.path.exists(log_file_path):
#         separator = "=" * 100
#
#         flog = open(log_file_path, 'w')
#         flog.write(separator + "\n")
#         flog.write("Read pair = " + read_pair_id + "\n")
#         flog.write(separator + "\n\n")
#
#         flog.flush()
#     else:
#         flog = open(log_file_path, 'a')
#
#     #-------------------------------------------------------------------------------------------------------------------
#     # Find contaminants in extended fragments
#     #-------------------------------------------------------------------------------------------------------------------
#
#     separator = "-" * 100
#
#     flash_file_name = read_pair_id + ".extendedFrags.fastq"
#     flash_file_path = os.path.join(flash_dir, flash_file_name)
#
#     sam_file_name = read_pair_id + "_" + db_name + "_Extended_Frags.sam"
#     sam_file_path = os.path.join(decon_dir, sam_file_name)
#
#     flog.write(separator + "\n")
#     flog.write("Contaminant DB = " + db_name + ", Extended Fragments\n")
#     flog.write(separator + "\n\n")
#
#     flog.flush()
#
#     the_args = [m_config['BW2_EXECUTABLE'], "-p", num_threads, "-x", contaminant_db_index_base, "-q", flash_file_path, "-S", sam_file_path, "--very-sensitive-local"]
#     the_cmd = " ".join(the_args)
#
#     subprocess.call(the_cmd, shell=True, stderr=flog)
#     flog.write("\n")
#
#     #-------------------------------------------------------------------------------------------------------------------
#     # Find contaminants in non-extended forward (R1) fragments
#     #-------------------------------------------------------------------------------------------------------------------
#
#     flash_file_name = read_pair_id + ".notCombined_1.fastq"
#     flash_file_path = os.path.join(flash_dir, flash_file_name)
#
#     sam_file_name = read_pair_id + "_" + db_name + "_Not_Combined_1.sam"
#     sam_file_path = os.path.join(decon_dir, sam_file_name)
#
#     flog.write(separator + "\n")
#     flog.write("Contaminant DB = " + db_name + ", Not Combined 1\n")
#     flog.write(separator + "\n\n")
#
#     flog.flush()
#
#     the_args = [m_config['BW2_EXECUTABLE'], "-p", num_threads, "-x", contaminant_db_index_base, "-q", flash_file_path, "-S", sam_file_path, "--very-sensitive-local"]
#     the_cmd = " ".join(the_args)
#
#     subprocess.call(the_cmd, shell=True, stderr=flog)
#     flog.write("\n")
#
#     #-------------------------------------------------------------------------------------------------------------------
#     # Find contaminants from non-extended reverse (R2) fragments
#     #-------------------------------------------------------------------------------------------------------------------
#
#     flash_file_name = read_pair_id + ".notCombined_2.fastq"
#     flash_file_path = os.path.join(flash_dir, flash_file_name)
#
#     sam_file_name = read_pair_id + "_" + db_name + "_Not_Combined_2.sam"
#     sam_file_path = os.path.join(decon_dir, sam_file_name)
#
#     flog.write(separator + "\n")
#     flog.write("Contaminant DB = " + db_name + ", Not Combined 2\n")
#     flog.write(separator + "\n\n")
#
#     flog.flush()
#
#     the_args = [m_config['BW2_EXECUTABLE'], "-p", num_threads, "-x", contaminant_db_index_base, "-q", flash_file_path, "-S", sam_file_path, "--very-sensitive-local"]
#     the_cmd = " ".join(the_args)
#
#     subprocess.call(the_cmd, shell=True, stderr=flog)
#     flog.write("\n")
#
#     #-------------------------------------------------------------------------------------------------------------------
#
#     flog.close()
#
#     end_time = time.time()
#     # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


def remove_contaminants(read_pair_id, db_name):
    start_time = time.time()
    decon_dir = get_decon_dir(read_pair_id)

    #-------------------------------------------------------------------------------------------------------------------
    # Create or open the  log file
    #-------------------------------------------------------------------------------------------------------------------

    log_file_name = read_pair_id + "_Decon_Contaminant_Removal.log"
    log_file_path = os.path.join(decon_dir, log_file_name)

    if not os.path.exists(log_file_path):
        separator = "=" * 100

        flog = open(log_file_path, 'w')
        flog.write(separator + "\n")
        flog.write("Read pair = " + read_pair_id + "\n")
        flog.write(separator + "\n\n")

        flog.flush()
    else:
        flog = open(log_file_path, 'a')

    #-------------------------------------------------------------------------------------------------------------------
    # Remove contaminants from extended fragments
    #-------------------------------------------------------------------------------------------------------------------

    separator = "-" * 100

    bw2_sam_file_name = read_pair_id + "_" + db_name + "_Extended_Frags.sam"
    bw2_sam_file_path = os.path.join(decon_dir, bw2_sam_file_name)

    picard_sam_file_name = read_pair_id + "_" + db_name + "_Unaligned_Extended_Frags.sam"
    picard_sam_file_path = os.path.join(decon_dir, picard_sam_file_name)

    flog.write(separator + "\n")
    flog.write("Contaminant DB = " + db_name + ", Extended Fragments\n")
    flog.write(separator + "\n\n")

    flog.flush()

    the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Unaligned", "I=" + bw2_sam_file_path, ">", picard_sam_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True,  stderr=flog)

    picard_fastq_file_name = read_pair_id + "_" + db_name + "_Unaligned_Extended_Frags.fastq"
    picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

    the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True,  stderr=flog)
    flog.write("\n")

    #-------------------------------------------------------------------------------------------------------------------
    # Remove contaminants from non-extended forward (R1) fragments
    #-------------------------------------------------------------------------------------------------------------------

    bw2_sam_file_name = read_pair_id + "_" + db_name + "_Not_Combined_1.sam"
    bw2_sam_file_path = os.path.join(decon_dir, bw2_sam_file_name)

    picard_sam_file_name = read_pair_id + "_" + db_name + "_Unaligned_Not_Combined_1.sam"
    picard_sam_file_path = os.path.join(decon_dir, picard_sam_file_name)

    flog.write(separator + "\n")
    flog.write("Contaminant DB = " + db_name + ", Not Combined 1\n")
    flog.write(separator + "\n\n")

    flog.flush()

    the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Unaligned", "I=" + bw2_sam_file_path, ">", picard_sam_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True,  stderr=flog)

    picard_fastq_file_name = read_pair_id + "_" + db_name + "_Unaligned_Not_Combined_1.fastq"
    picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

    the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True,  stderr=flog)
    flog.write("\n")

    #-------------------------------------------------------------------------------------------------------------------
    # Remove contaminants from non-extended reverse (R2) fragments
    #-------------------------------------------------------------------------------------------------------------------

    bw2_sam_file_name = read_pair_id + "_" + db_name + "_Not_Combined_2.sam"
    bw2_sam_file_path = os.path.join(decon_dir, bw2_sam_file_name)

    picard_sam_file_name = read_pair_id + "_" + db_name + "_Unaligned_Not_Combined_2.sam"
    picard_sam_file_path = os.path.join(decon_dir, picard_sam_file_name)

    flog.write(separator + "\n")
    flog.write("Contaminant DB = " + db_name + ", Not Combined 2\n")
    flog.write(separator + "\n\n")

    flog.flush()

    the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Unaligned", "I=" + bw2_sam_file_path, ">", picard_sam_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True,  stderr=flog)

    picard_fastq_file_name = read_pair_id + "_" + db_name + "_Unaligned_Not_Combined_2.fastq"
    picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

    the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True,  stderr=flog)
    flog.write("\n")

    flog.close()

    end_time = time.time()

    # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


def concat_decon_files(read_pair_id):
    decon_dir = get_decon_dir(read_pair_id)

    fastq_base_file_names = ["Unaligned_Extended_Frags", "Unaligned_Not_Combined_1", "Unaligned_Not_Combined_2"]
    for base_file_name in fastq_base_file_names:
        fastq_file_pattern = read_pair_id + "_*_" + base_file_name + ".fastq"
        fastq_file_pattern = os.path.join(decon_dir, fastq_file_pattern)
        fastq_files = glob.glob(fastq_file_pattern)

        cat_fastq_file_name = read_pair_id + "_" + base_file_name + "_Concat.fastq"
        cat_fastq_file_path = os.path.join(decon_dir, cat_fastq_file_name)

        with open(cat_fastq_file_path, 'w') as fout:
            for fname in fastq_files:
                with open(fname) as fin:
                    for line in fin:
                        fout.write(line)


def decon(read_pair_id):
    #-------------------------------------------------------------------------------------------------------------------
    # Create the decon output directory if it doesn't exist
    #-------------------------------------------------------------------------------------------------------------------

    decon_dir = get_decon_dir(read_pair_id)

    if os.path.exists(decon_dir):
        shutil.rmtree(decon_dir)

    os.makedirs(decon_dir)

    #-------------------------------------------------------------------------------------------------------------------
    # Check if the contaminant databases exist and build them if they don't
    #-------------------------------------------------------------------------------------------------------------------

    contaminant_dbs = m_param['decon']['contaminant_dbs']
    contaminant_dbs = contaminant_dbs.split(",")

    flash_dir = get_flash_dir(read_pair_id)

    step_id = 1
    final_step = False

    for db_name in contaminant_dbs:
        db_name = db_name.strip()
        db_exists = contaminant_db_exists(db_name)
        if not db_exists:
            build_contaminant_db(db_name)

        if step_id == len(contaminant_dbs):
            final_step = True

        find_contaminants(read_pair_id, db_name, step_id, final_step)
        # remove_contaminants(read_pair_id, db_name)

        step_id += 1

    # concat_decon_files(read_pair_id)


def trim_reads(read_pair_id):
   # print("Removing trimming reads using Trimmomatic..")

    start_time = time.time()
    decon_dir = get_decon_dir(read_pair_id)
    trimmomatic_dir = get_trimmomatic_dir(read_pair_id)

    if os.path.exists(trimmomatic_dir):
        shutil.rmtree(trimmomatic_dir)

    os.makedirs(trimmomatic_dir)

    contaminant_dbs = m_param['decon']['contaminant_dbs']
    contaminant_dbs = contaminant_dbs.split(",")

   #-------------------------------------------------------------------------------------------------------------------
    # Create or open the  log file
    #-------------------------------------------------------------------------------------------------------------------

    log_file_name = read_pair_id + "_Trimmomatic.log"
    log_file_path = os.path.join(trimmomatic_dir, log_file_name)

    if not os.path.exists(log_file_path):
        separator = "=" * 100

        flog = open(log_file_path, 'w')
        flog.write(separator + "\n")
        flog.write("Read pair = " + read_pair_id + "\n")
        flog.write(separator + "\n\n")

        flog.flush()
    else:
        flog = open(log_file_path, 'a')

    #-------------------------------------------------------------------------------------------------------------------
    # Trimming extended fragments
    #-------------------------------------------------------------------------------------------------------------------

    for db_name in contaminant_dbs:
        the_file_pattern = "Final_*_Unaligned_Extended_Frags.fastq"
        the_file_pattern = os.path.join(decon_dir, the_file_pattern)
        picard_fastq_file_path = glob.glob(the_file_pattern)
        picard_fastq_file_path = picard_fastq_file_path[0]

        # trimmomatic_adapters_path = "/home/whit040/Programs/Trimmomatic-0.33/adapters/TruSeq2-SE.fa"
        trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.fastq"
        trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

        trimmomatic_log_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.log"
        trimmomatic_log_file_path = os.path.join(trimmomatic_dir, trimmomatic_log_file_name)

        flog.write(separator + "\n")
        flog.write("Contaminant DB = " + db_name + ", Extended Fragments\n")
        flog.write(separator + "\n\n")

        flog.flush()

        the_args = [m_config['TRIMMOMATIC_EXECUTABLE'], "SE", "-phred33", "-trimlog", trimmomatic_log_file_path, picard_fastq_file_path, trimmomatic_fastq_file_path,
                    "ILLUMINACLIP:" + m_config['TRIMMOMATIC_ADAPTER_FILE'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:50",
                    ""]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True, stderr=flog)

        flog.write("\n")

        #-------------------------------------------------------------------------------------------------------------------
        # Trimming non-extended forward (R1) fragments
        #-------------------------------------------------------------------------------------------------------------------

        the_file_pattern = "Final_*_Unaligned_Not_Combined_1.fastq"
        the_file_pattern = os.path.join(decon_dir, the_file_pattern)
        picard_fastq_file_path = glob.glob(the_file_pattern)
        picard_fastq_file_path = picard_fastq_file_path[0]

        # picard_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Concat.fastq"
        # picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

        trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.fastq"
        trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

        trimmomatic_log_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.log"
        trimmomatic_log_file_path = os.path.join(trimmomatic_dir, trimmomatic_log_file_name)

        flog.write(separator + "\n")
        flog.write("Contaminant DB = " + db_name + ", Not Combined 1\n")
        flog.write(separator + "\n\n")

        flog.flush()

        the_args = [m_config['TRIMMOMATIC_EXECUTABLE'], "SE", "-phred33", "-trimlog", trimmomatic_log_file_path, picard_fastq_file_path, trimmomatic_fastq_file_path,
                    "ILLUMINACLIP:" + m_config['TRIMMOMATIC_ADAPTER_FILE'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:50"]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True, stderr=flog)

        flog.write("\n")

        #-------------------------------------------------------------------------------------------------------------------
        # Trimming non-extended reverse (R2) fragments
        #-------------------------------------------------------------------------------------------------------------------

        # picard_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_2_Concat.fastq"
        # picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

        the_file_pattern = "Final_*_Unaligned_Not_Combined_2.fastq"
        the_file_pattern = os.path.join(decon_dir, the_file_pattern)
        picard_fastq_file_path = glob.glob(the_file_pattern)
        picard_fastq_file_path = picard_fastq_file_path[0]

        trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_2_Trimmed.fastq"
        trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

        trimmomatic_log_file_name = read_pair_id + "_Unaligned_Not_Combined_2_Trimmed.log"
        trimmomatic_log_file_path = os.path.join(trimmomatic_dir, trimmomatic_log_file_name)

        flog.write(separator + "\n")
        flog.write("Contaminant DB = " + db_name + ", Not Combined 2\n")
        flog.write(separator + "\n\n")

        flog.flush()

        the_args = [m_config['TRIMMOMATIC_EXECUTABLE'], "SE", "-phred33", "-trimlog", trimmomatic_log_file_path, picard_fastq_file_path, trimmomatic_fastq_file_path,
                    "ILLUMINACLIP:" + m_config['TRIMMOMATIC_ADAPTER_FILE'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:50"]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True, stderr=flog)

        flog.write("\n")

    flog.close()

    end_time = time.time()
    # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


def qc_reads(read_pair_id):
    # print("Generating Quality control statistics using FastQC..")

    start_time = time.time()

    trimmomatic_dir = get_trimmomatic_dir(read_pair_id)
    fastqc_dir = get_fastqc_dir(read_pair_id)

    if os.path.exists(fastqc_dir):
        shutil.rmtree(fastqc_dir)

    os.makedirs(fastqc_dir)

    flog = open(os.devnull, 'w')

    #-------------------------------------------------------------------------------------------------------------------
    # QC statistics for extended fragments
    #-------------------------------------------------------------------------------------------------------------------

    trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.fastq"
    trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

    the_args = [m_config['FASTQC_EXECUTABLE'], trimmomatic_fastq_file_path, "-o", fastqc_dir]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)

    #-------------------------------------------------------------------------------------------------------------------
    # QC statistics for non-extended forward (R1) fragments
    #-------------------------------------------------------------------------------------------------------------------

    trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.fastq"
    trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

    the_args = [m_config['FASTQC_EXECUTABLE'], trimmomatic_fastq_file_path, "-o", fastqc_dir]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)

    #-------------------------------------------------------------------------------------------------------------------
    # QC statistics for non-extended reverse (R2) fragments
    #-------------------------------------------------------------------------------------------------------------------

    trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_2_Trimmed.fastq"
    trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

    the_args = [m_config['FASTQC_EXECUTABLE'], trimmomatic_fastq_file_path, "-o", fastqc_dir]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)

    flog.close()

    end_time = time.time()

    # print("\tCompleted in " + str(end_time - start_time) + " seconds.")


def interleave_reads(read_pair_id):
    # print("Generating Quality control statistics using FastQC..")
    trimmomatic_dir = get_trimmomatic_dir(read_pair_id)
    interleave_dir = get_interleave_dir(read_pair_id)

    interleaved_file_name = read_pair_id + "_Trimmed_Interleaved.fastq"
    interleaved_file_path = os.path.join(interleave_dir, interleaved_file_name)

    if os.path.exists(interleave_dir):
        shutil.rmtree(interleave_dir)

    os.makedirs(interleave_dir)

    r1_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.fastq"
    r1_file_path = os.path.join(trimmomatic_dir, r1_file_name)

    r2_file_name = read_pair_id + "_Unaligned_Not_Combined_2_Trimmed.fastq"
    r2_file_path = os.path.join(trimmomatic_dir, r2_file_name)

    fr1 = open(r1_file_path, 'r')
    fr2 = open(r2_file_path, 'r')
    fout = open(interleaved_file_path, 'w')

    while True:
        r1_id = fr1.readline()
        if not r1_id: break
        r1_seq = fr1.readline()
        r1_plus = fr1.readline()
        r1_quals = fr1.readline()

        r2_id = fr2.readline()
        r2_seq = fr2.readline()
        r2_plus = fr2.readline()
        r2_quals = fr2.readline()

        fout.write(r1_id)
        fout.write(r1_seq)
        fout.write(r1_plus)
        fout.write(r1_quals)

        fout.write(r2_id)
        fout.write(r2_seq)
        fout.write(r2_plus)
        fout.write(r2_quals)

    fr1.close()
    fr2.close()
    fout.close()


def get_forward_read_path(read_pair_id):
    file_name = read_pair_id + "_R1.fastq"
    file_path = os.path.join(m_config['INPUT_DIR'], file_name)

    return file_path


def get_reverse_read_path(read_pair_id):
    file_name = read_pair_id + "_R2.fastq"
    file_path = os.path.join(m_config['INPUT_DIR'], file_name)

    return file_path


def get_decon_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Decon")

    return the_dir


def get_flash_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Flash")

    return the_dir


def get_trimmomatic_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Trimmomatic")

    return the_dir


def get_fastqc_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "FastQC")

    return the_dir


def get_interleave_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Interleave")

    return the_dir


# def create_adapter_file():
#     input_file_path = "/home/whit040/Desktop/contaminant_list.txt"
#     output_file_path = "/home/whit040/Programs/Trimmomatic-0.33/adapters/Adapters.fa"
#
#     fin = open(input_file_path, 'r')
#     fout = open(output_file_path, 'w')
#
#     for line in fin:
#         line = line.strip()
#         if re.match("#", line) or line == "":
#             continue
#
#         line_parts = line.split("\t")
#         fout.write(">" + line_parts[0] + "\n")
#
#         for p in line_parts[1:]:
#             p = p.strip()
#             if p != "":
#                 fout.write(p + "\n")
#                 break
#
#     fout.close()
#     fin.close()


def worker(args):
    read_pair_id = args

    merge_reads(read_pair_id)
    decon(read_pair_id)
    trim_reads(read_pair_id)
    qc_reads(read_pair_id)
    interleave_reads(read_pair_id)


def run_parallel():
    #---------------------------------------------------------------------------------------------------------------
    # Build the contaminant databases if they don't exist
    #---------------------------------------------------------------------------------------------------------------

    contaminant_dbs = m_param['decon']['contaminant_dbs']
    contaminant_dbs = contaminant_dbs.split(",")

    for db_name in contaminant_dbs:
        db_name = db_name.strip()
        db_exists = contaminant_db_exists(db_name)
        if not db_exists:
            build_contaminant_db(db_name)

    #---------------------------------------------------------------------------------------------------------------
    # Create the job queue and then run it
    #---------------------------------------------------------------------------------------------------------------

    num_procs = 5

    if num_procs > multiprocessing.cpu_count() - 1:
        num_procs = multiprocessing.cpu_count() - 1

    task_args = []
    for read_pair_id in m_read_pairs.keys():
        read_pair_output_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
        if os.path.exists(read_pair_output_dir):
            print("Skipping read pair '" + read_pair_id + "'. Already exists.")
            continue

        task_args.append(read_pair_id,)

    p = multiprocessing.Pool(processes=num_procs)
    p.map(worker, task_args)
    p.close()
    p.join()


def run_serial():
    performance_log = {}

    seperator = "-" * 100

    for read_pair_id in m_read_pairs.keys():
        performance_log[read_pair_id] = {}

        print(seperator)
        print("Working on read pair '" + read_pair_id + "'")
        print(seperator + "\n")

        #---------------------------------------------------------------------------------------------------------------
        # Merge reads
        #---------------------------------------------------------------------------------------------------------------

        print("Merging reads...")
        start_time = time.time()
        merge_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['merging'] = end_time - start_time

        #---------------------------------------------------------------------------------------------------------------
        # Decon reads
        #---------------------------------------------------------------------------------------------------------------

        print("Decontaminating reads...")
        start_time = time.time()
        decon(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['decon'] = end_time - start_time

        #---------------------------------------------------------------------------------------------------------------
        # Trim reads
        #---------------------------------------------------------------------------------------------------------------

        print("Trimming reads...")
        start_time = time.time()
        trim_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['trimming'] = end_time - start_time

        #---------------------------------------------------------------------------------------------------------------
        # QC reads
        #---------------------------------------------------------------------------------------------------------------

        print("QC reads...")
        start_time = time.time()
        qc_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['qc'] = end_time - start_time

        #---------------------------------------------------------------------------------------------------------------
        # Interleave reads
        #---------------------------------------------------------------------------------------------------------------

        print("Interleave reads...")
        start_time = time.time()
        interleave_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['interleaving'] = end_time - start_time
        print("\n")

    print(seperator)
    print("Performance")
    print(seperator + "\n\n")

    print(performance_log)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_dir', help='The directory containing the sequence files (read pairs)')
    parser.add_argument('-output_dir', help='The directory to use for writing out the results.')
    parser.add_argument('-config_file', help='The path to the config file.')
    parser.add_argument('-param_file', help='The path to the parameter file.')

    args = parser.parse_args()

    global m_input_dir
    m_input_dir = args.input_dir

    global m_output_dir
    m_output_dir = args.output_dir

    global m_config_file
    m_config_file = args.config_file

    global m_param_file
    m_param_file = args.param_file

    init_run_info()

    # run_serial()
    run_parallel()


if __name__ == "__main__":
    main()
