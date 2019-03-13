import argparse
import glob
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from psutil import virtual_memory

from scripts.Chunker import Chunker
from scripts import function_parser

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
        if program_name not in m_param.keys():
            m_param[program_name] = {}

        line_data = line_data[1].split(" ")
        param_name = line_data[0]
        param_val = line_data[1:]
        param_val = " ".join(param_val)

        m_param[program_name][param_name] = param_val

    fIn.close()


def build_dbs(database_dir):
    # keys = ['*.bck', '*.des', '*-names.txt', '*.prj', '*.sds', '*.ssp', '*.suf', '*.tis']

    dbs = []
    if os.path.exists(database_dir):
        dbs.extend(glob.glob(os.path.join(database_dir, '*')))
        dbs = [x for x in dbs if os.path.isfile(x)]
    else:
        raise("Error: no %s database folder." % os.path.basename(database_dir))

    formatted_dir = os.path.join(database_dir, 'formatted')
    if not os.path.exists(formatted_dir):
        os.makedirs(formatted_dir)

    if len(os.listdir(formatted_dir)) == 0:
        return

    for db in dbs:
        db_name = os.path.basename(os.path.splitext(db)[0] + '.dmnd')
        # lastdb+ call
        if os.path.basename(database_dir) == 'functional':
            cmd = '%s --makedb --in %s -d %s' % (m_config['DIAMOND_EXECUTABLE'], db, os.path.join(formatted_dir, db_name))
        else:
            cmd = '%s %s %s' % (m_config['LASTDB+_EXECUTABLE'], db, db)
        with open(os.devnull, 'w') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=f)


def build_annotation_dbs():
    functional = os.path.join(m_config['ANNOTATION_DB_DIR'], 'functional')
    # taxonomic = os.path.join(m_config['ANNOTATION_DB_DIR'], 'taxonomic')

    build_dbs(functional)
    # build_dbs(taxonomic)


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

    #
    # Create or open the  log file
    #
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

    #
    # Find contaminants in extended fragments
    #
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

        picard_sam_file_name = step_prefix + read_pair_id + "_" + db_name + "_Unaligned_" + base_names[2] + ".sam"
        picard_sam_file_path = os.path.join(decon_dir, picard_sam_file_name)

        the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Unaligned", "I=" + bw2_output_file_path, ">", picard_sam_file_path]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True, stderr=flog)

        picard_fastq_file_name = step_prefix + read_pair_id + "_" + db_name + "_Unaligned_" + base_names[2] + ".fastq"
        picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

        the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
        the_cmd = " ".join(the_args)
        subprocess.call(the_cmd, shell=True, stderr=flog)

        flog.write("\n")


def remove_contaminants(read_pair_id, db_name):
    start_time = time.time()
    decon_dir = get_decon_dir(read_pair_id)

    #
    # Create or open the  log file
    #
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

    #
    # Remove contaminants from extended fragments
    #
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
    subprocess.call(the_cmd, shell=True, stderr=flog)

    picard_fastq_file_name = read_pair_id + "_" + db_name + "_Unaligned_Extended_Frags.fastq"
    picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

    the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stderr=flog)
    flog.write("\n")

    #
    # Remove contaminants from non-extended forward (R1) fragments
    #
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
    subprocess.call(the_cmd, shell=True, stderr=flog)

    picard_fastq_file_name = read_pair_id + "_" + db_name + "_Unaligned_Not_Combined_1.fastq"
    picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

    the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stderr=flog)
    flog.write("\n")

    #
    # Remove contaminants from non-extended reverse (R2) fragments
    #
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
    subprocess.call(the_cmd, shell=True, stderr=flog)

    picard_fastq_file_name = read_pair_id + "_" + db_name + "_Unaligned_Not_Combined_2.fastq"
    picard_fastq_file_path = os.path.join(decon_dir, picard_fastq_file_name)

    the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + picard_sam_file_path, "F=" + picard_fastq_file_path]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stderr=flog)
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
    #
    # Create the decon output directory if it doesn't exist
    #
    decon_dir = get_decon_dir(read_pair_id)

    if os.path.exists(decon_dir):
        shutil.rmtree(decon_dir)

    os.makedirs(decon_dir)

    #
    # Check if the contaminant databases exist and build them if they don't
    #
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

    #
    # Create or open the  log file
    #
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

    #
    # Trimming extended fragments
    #
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

        #
        # Trimming non-extended forward (R1) fragments
        #
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

        #
        # Trimming non-extended reverse (R2) fragments
        #

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

    #
    # QC statistics for extended fragments
    #
    trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.fastq"
    trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

    the_args = [m_config['FASTQC_EXECUTABLE'], trimmomatic_fastq_file_path, "-o", fastqc_dir]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)

    #
    # QC statistics for non-extended forward (R1) fragments
    #
    trimmomatic_fastq_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.fastq"
    trimmomatic_fastq_file_path = os.path.join(trimmomatic_dir, trimmomatic_fastq_file_name)

    the_args = [m_config['FASTQC_EXECUTABLE'], trimmomatic_fastq_file_path, "-o", fastqc_dir]
    the_cmd = " ".join(the_args)
    subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)

    #
    # QC statistics for non-extended reverse (R2) fragments
    #
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
        if not r1_id:
            break
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


def merge_trim_outputs(read_pair_id):
    # relevant input directories
    trimmomatic_dir = get_trimmomatic_dir(read_pair_id)
    interleave_dir = get_interleave_dir(read_pair_id)
    assembly_dir = get_assembly_dir(read_pair_id)

    # trimmomatic file
    trimmomatic_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.fastq"
    trimmomatic_file_path = os.path.join(trimmomatic_dir, trimmomatic_file_name)

    # interleave file
    interleaved_file_name = read_pair_id + "_Trimmed_Interleaved.fastq"
    interleaved_file_path = os.path.join(interleave_dir, interleaved_file_name)

    # assembly file (new file!)
    assembly_file_name = read_pair_id + "_Ext-IL_Trimmed.fastq"
    assembly_file_path = os.path.join(assembly_dir, assembly_file_name)

    # create assembly folder
    if os.path.exists(assembly_dir):
        shutil.rmtree(assembly_dir)

    os.makedirs(assembly_dir)

    # cat trimmomatic and interleave into new file
    with open(assembly_file_path, 'w') as outfile:
        for fname in [trimmomatic_file_path, interleaved_file_path]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def megahit(read_pair_id):
    assembly_dir = get_assembly_dir(read_pair_id)

    assembly_file_name = read_pair_id + "_Ext-IL_Trimmed.fastq"
    assembly_file_path = os.path.join(assembly_dir, assembly_file_name)

    # megahit dir
    megahit_dir = os.path.join(assembly_dir, 'MegaHit')

    # load merged trim file into megahit
    the_cmd = '%s -m %s -l %s -r %s --k-min %s --k-max %s --out-dir %s' % \
              (m_config['MEGAHIT_EXECUTABLE'], m_param['megahit']['max_mem'], m_param['megahit']['length_of_library_insert'],
               assembly_file_path, m_param['megahit']['kmer_min'], m_param['megahit']['kmer_max'], megahit_dir)

    with open(os.path.join(assembly_dir, 'megahit_log.txt'), 'w') as flog:
        subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)

    # move output file to assembly root
    # might be "contigs" below
    shutil.copy2(os.path.join(megahit_dir, 'final.contigs.fa'), os.path.join(assembly_dir, read_pair_id + '_final_contigs.fasta'))


def subsample(read_pair_id):
    # subsampling_length.py implementation to get >1K reads (user specified length)
    assembly_dir = get_assembly_dir(read_pair_id)

    assembly_name = read_pair_id + '_final_contigs.fasta'
    assembly_path = os.path.join(assembly_dir, assembly_name)

    outfile_name = assembly_name[:-6] + '_1k.fasta'
    outfile_path = os.path.join(assembly_dir, outfile_name)

    the_cmd = 'python ./scripts/subsampler_length.py -f %s -m %s -n %s > %s' % \
              (assembly_path, m_param['subsampler']['min'], m_param['subsampler']['max'], outfile_path)

    with open(os.devnull, 'w') as flog:
        subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)


def build_assembly_db(read_pair_id):
    assembly_dir = get_assembly_dir(read_pair_id)
    coverage_dir = get_coverage_dir(read_pair_id)

    if os.path.exists(coverage_dir):
        shutil.rmtree(coverage_dir)

    os.makedirs(coverage_dir)

    coverage_dir = get_coverage_dir(read_pair_id)

    assembly_name = read_pair_id + '_final_contigs.fasta'
    assembly_path = os.path.join(assembly_dir, assembly_name)

    subsampled_name = assembly_name[:-6] + '_1k.fasta'
    subsampled_path = os.path.join(assembly_dir, subsampled_name)

    db_name = os.path.splitext(subsampled_name)[0]
    db_path = os.path.join(coverage_dir, db_name)

    the_args = [m_config['BW2_BUILD_EXECUTABLE'], subsampled_path, db_path]
    # the_cmd = " ".join(the_args)

    flog = os.path.join(coverage_dir, 'assembly_db_format.log')
    with open(flog, 'w') as f:
        subprocess.call(the_args, shell=False, stderr=f, stdout=f)


def map_to_assembly(read_pair_id):
    trimmomatic_dir = get_trimmomatic_dir(read_pair_id)
    assembly_dir = get_assembly_dir(read_pair_id)
    coverage_dir = get_coverage_dir(read_pair_id)

    # get db base name
    assembly_name = read_pair_id + '_final_contigs.fasta'

    assembly_path = os.path.join(assembly_dir, assembly_name)

    db_name = assembly_name[:-6] + '_1k'
    db_path = os.path.join(coverage_dir, db_name)

    # get 5 read files
    r1_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.fastq"
    r1_file_path = os.path.join(trimmomatic_dir, r1_file_name)

    r2_file_name = read_pair_id + "_Unaligned_Not_Combined_2_Trimmed.fastq"
    r2_file_path = os.path.join(trimmomatic_dir, r2_file_name)

    extended_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.fastq"
    extended_file_path = os.path.join(trimmomatic_dir, extended_file_name)

    ext_il_name = read_pair_id + "_Ext-IL_Trimmed.fastq"
    ext_il_path = os.path.join(assembly_dir, ext_il_name)

    # cat R1 and extended
    extended_cat_r1_name = read_pair_id + "_Unaligned_Ext_R1_Trimmed.fastq"
    extended_cat_r1_path = os.path.join(coverage_dir, extended_cat_r1_name)

    with open(extended_cat_r1_path, 'w') as outfile:
        for fname in [extended_file_path, r1_file_path]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    # call bw2 for all 5 files
    jobs = []
    jobs_picard = []
    num_threads = multiprocessing.cpu_count()
    with open(os.path.join(coverage_dir, 'bw2.log'), 'w') as flog:
        for inF in [r1_file_path, r2_file_path, extended_file_path, extended_cat_r1_path, ext_il_path]:

            bn = os.path.splitext(os.path.basename(inF))[0]
            outF = os.path.join(coverage_dir, bn + '.sam')

            the_args = [m_config['BW2_EXECUTABLE'], "-p", str(num_threads), "-x", db_path, "-q", inF, "-S", outF, "--very-sensitive-local"]
            the_args = " ".join(the_args)
            p = subprocess.Popen(the_args, shell=True, stderr=flog)
            jobs.append(p)
            jobs_picard.append(outF)

    # wait for all bw2 calls to finish
    for p in jobs:
        p.wait()

    with open(os.path.join(coverage_dir, 'picard.log'), 'w') as flog:
        for inF in jobs_picard:

            outF_ua = os.path.splitext(inF)[0] + '_Unaligned'
            outF_a = os.path.splitext(inF)[0] + '_Aligned'

            # run picard with unaligned
            the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Unaligned", "I=" + inF, ">", outF_ua + '.sam']
            the_cmd = " ".join(the_args)
            subprocess.call(the_cmd, shell=True, stderr=flog, stdout=flog)

            # run picard with aligned
            the_args = [m_config['PICARD_VIEW_SAM_EXECUTABLE'], "ALIGNMENT_STATUS=Aligned", "I=" + inF, ">", outF_a + '.sam']
            the_cmd = " ".join(the_args)
            subprocess.call(the_cmd, shell=True, stderr=flog, stdout=flog)

            # convert unaligned .sam files to .bam files
            the_args = [m_config['PICARD_SAM_TO_FASTQ_EXECUTABLE'], "I=" + outF_ua + '.sam', "F=" + outF_ua + '.fastq']
            the_cmd = " ".join(the_args)
            subprocess.call(the_cmd, shell=True, stderr=flog, stdout=flog)

            # convert aligned .sam files to .bam files
            the_args = [m_config['SAMTOOLS_EXECUTABLE'], 'view', '-bS', outF_a + '.sam', ">", outF_a + '.bam']
            the_cmd = " ".join(the_args)
            subprocess.call(the_cmd, shell=True, stderr=flog, stdout=flog)


def get_stats(read_pair_id):
    assembly_dir = get_assembly_dir(read_pair_id)
    stats_dir = get_stats_dir(read_pair_id)

    assembly_name = read_pair_id + '_final_contigs.fasta'
    assembly_path = os.path.join(assembly_dir, assembly_name)

    # create binning folder
    if os.path.exists(stats_dir):
        shutil.rmtree(stats_dir)

    os.mkdir(stats_dir)

    # find subsampled file
    subsampled_file_name = assembly_name[:-6] + '_1k.fasta'
    subsampled_file_path = os.path.join(assembly_dir, subsampled_file_name)

    subsampled_stats_file_name = subsampled_file_name[:-6] + '_stats.txt'
    subsampled_stats_file_path = os.path.join(stats_dir, subsampled_stats_file_name)

    final_contigs_stats_file_name = assembly_name[:-6] + '_stats.txt'
    final_contigs_stats_file_path = os.path.join(stats_dir, final_contigs_stats_file_name)

    # run it on ouput of subsample
    subsampled_perl_cmd = 'perl ./scripts/CountFasta.pl %s > %s' % (subsampled_file_path, subsampled_stats_file_path)
    final_contigs_perl_cmd = 'perl ./scripts/CountFasta.pl %s > %s' % (assembly_path, final_contigs_stats_file_path)

    with open(os.path.join(stats_dir, 'stats.log'), 'w') as flog:
        subprocess.call(subsampled_perl_cmd, shell=True, stdout=flog, stderr=flog)
        subprocess.call(final_contigs_perl_cmd, shell=True, stdout=flog, stderr=flog)


def maxbin(read_pair_id):
    # get relevant directories
    assembly_dir = get_assembly_dir(read_pair_id)
    trimmomatic_dir = get_trimmomatic_dir(read_pair_id)
    binning_dir = get_binning_dir(read_pair_id)

    # create binning folder
    if os.path.exists(binning_dir):
        shutil.rmtree(binning_dir)

    os.mkdir(binning_dir)

    # trimmomatic files
    trimmomatic_file_name = read_pair_id + "_Unaligned_Extended_Frags_Trimmed.fastq"
    trimmomatic_file_path = os.path.join(trimmomatic_dir, trimmomatic_file_name)

    trimmomatic2_file_name = read_pair_id + "_Unaligned_Not_Combined_1_Trimmed.fastq"
    trimmomatic2_file_path = os.path.join(trimmomatic_dir, trimmomatic2_file_name)

    # output of cat
    catted_trims_name = read_pair_id + '_Extended_Frags_Trimmed_R1.fastq'
    catted_trims_path = os.path.join(binning_dir, catted_trims_name)

    # cat trimmomatic and interleave into new file
    with open(catted_trims_path, 'w') as outfile:
        for fname in [trimmomatic_file_path, trimmomatic2_file_path]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    assembly_name = read_pair_id + '_final_contigs.fasta'
    assembly_path = os.path.join(assembly_dir, assembly_name)

    subsampled_file_name = assembly_name[:-6] + '_1k.fasta'
    subsampled_file_path = os.path.join(assembly_dir, subsampled_file_name)

    the_cmd = '%s -contig %s -reads %s -out %s -thread %s' % \
              (m_config['MAXBIN_EXECUTABLE'], subsampled_file_path, catted_trims_path, binning_dir, multiprocessing.cpu_count())

    with open(os.path.join(binning_dir, 'maxbin_log.txt'), 'w') as flog:
        subprocess.call(the_cmd, shell=True, stdout=flog, stderr=flog)


def prodigal(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    assembly_dir = get_assembly_dir(read_pair_id)

    # create annotation folder
    if os.path.exists(annotation_dir):
        shutil.rmtree(annotation_dir)

    os.mkdir(annotation_dir)
    orf_dir = os.path.join(annotation_dir, 'ORFprediction')

    if not os.path.exists(orf_dir):
        os.mkdir(orf_dir)

    assembly_name = read_pair_id + '_final_contigs_1k.fasta'
    assembly_path = os.path.join(assembly_dir, assembly_name)

    query = os.path.splitext(assembly_name)[0]

    gff_name = query + '.gff'
    gff_path = os.path.join(orf_dir, gff_name)

    proteins_name = query + '.faa'
    proteins_path = os.path.join(orf_dir, proteins_name)

    nuc_name = query + '.ffn'
    nuc_path = os.path.join(orf_dir, nuc_name)

    with open(os.devnull) as flog:
        cmd = 'prodigal -i %s -o %s -a %s -d %s -f gff -p meta' % (assembly_path, gff_path, proteins_path, nuc_path)
        subprocess.call(cmd, shell=True, stderr=flog, stdout=flog)


def diamond(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    assembly_dir = get_assembly_dir(read_pair_id)
    func_dir = os.path.join(annotation_dir, 'functional')

    if not os.path.exists(func_dir):
        os.mkdir(func_dir)

    orf_dir = os.path.join(annotation_dir, 'ORFprediction')

    assembly_name = read_pair_id + '_final_contigs_1k.fasta'
    assembly_path = os.path.join(assembly_dir, assembly_name)

    query = os.path.splitext(assembly_name)[0]
    proteins_name = query + '.faa'
    proteins_path = os.path.join(orf_dir, proteins_name)

    refdbs = glob.glob(os.path.join(m_config['ANNOTATION_DB_DIR'], 'functional', 'formatted', '*'))
    refdbs = [x for x in refdbs if os.path.isfile(x)]

    split_dir = os.path.join(func_dir, 'tmp/split/')
    out_dir = os.path.join(func_dir, 'tmp/out/')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    c = Chunker(proteins_path, split_dir, chunksize='5M', delim='>')

    for refdb in refdbs:
        outFiles = []
        db_name = os.path.splitext(os.path.basename(refdb))[0]
        with open(os.path.join(func_dir, 'diamond_%s.log' % db_name), 'w') as flog:
            for inF in c.files:
                outName = os.path.splitext(os.path.basename(inF))[0] + '_' + db_name
                outF = os.path.join(out_dir, outName)
                outFiles.append(outF)

                cmd = "diamond blastp --db %s --query %s --more-sensitive -f 6 -k 1 -b 1 -o %s" % (refdb, inF, outF)
                subprocess.call(cmd, shell=True, stderr=flog, stdout=flog)

        # cat diamond output together
        final_output = os.path.join(func_dir, "%s_diamond_%s" % (query, os.path.basename(os.path.splitext(refdb)[0])))

        cat_cmd = 'cat'
        for o in outFiles:
            cat_cmd += ' ' + o
        cat_cmd += ' > ' + final_output

        subprocess.call(cat_cmd, shell=True)

        # remove temp folders
        # shutil.rmtree(os.path.join(func_dir, 'tmp/'))


def sword(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    assembly_dir = get_assembly_dir(read_pair_id)

    tax_dir = os.path.join(annotation_dir, 'taxonomic')
    if not os.path.exists(tax_dir):
        os.mkdir(tax_dir)

    subsampled_name = read_pair_id + '_final_contigs_1k.fasta'
    subsampled_path = os.path.join(assembly_dir, subsampled_name)

    dbs = glob.glob(os.path.join(m_config['ANNOTATION_DB_DIR'], 'taxonomic', '*'))
    dbs = [x for x in dbs if os.path.isfile(x)]

    # params
    e_value = m_param['sword']['e_value']
    cores = m_param['sword']['num_threads']
    kmer_length = m_param['sword']['kmer_length']

    split_dir = os.path.join(tax_dir, 'tmp/split/')
    out_dir = os.path.join(tax_dir, 'tmp/out/')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    c = Chunker(subsampled_path, split_dir, chunksize='5M', delim='>')

    for db in dbs:
        outFiles = []
        db_name = os.path.splitext(os.path.basename(db))[0]
        with open(os.path.join(tax_dir, 'sword_%s.log' % db_name), 'w') as flog:
            for inF in c.files:
                outName = os.path.splitext(os.path.basename(inF))[0]
                outF = os.path.join(out_dir, outName)
                outFiles.append(outF)

                # construct command
                cmd = ['sword', '-i', inF, '-j', db, '-v', e_value, '-t', cores, '-k', kmer_length, '-o', outF]
                subprocess.call(cmd, stderr=flog, stdout=flog)

        # cat sword output together
        final_output = os.path.join(tax_dir, "%s_sword_%s" % (os.path.splitext(subsampled_name)[0], db_name))

        cat_cmd = 'cat'
        for o in outFiles:
            cat_cmd += ' ' + o
        cat_cmd += ' > ' + final_output
        subprocess.call(cat_cmd, shell=True)

        # remove temp folders
        # shutil.rmtree(os.path.join(tax_dir, 'tmp/'))


def blast_parsing(read_pair_id):
    # diamond output
    annotation_dir = get_annotation_dir(read_pair_id)
    func_dir = os.path.join(annotation_dir, 'functional')
    blast_dir = os.path.join(func_dir, 'blast_out')

    if not os.path.exists(blast_dir):
        os.mkdir(blast_dir)

    dbs = glob.glob(os.path.join(func_dir, '*'))
    dbs = [x for x in dbs if os.path.isfile(x)]
    dbs = [x for x in dbs if '.log' not in x]

    for db in dbs:
        db_sorted = os.path.splitext(db)[0] + '_sorted.tsv'
        db_name = os.path.basename(db_sorted)

        sort_cmd = ['sort', '-k1,1', '-k12,12rn', db, '>', db_sorted]
        sort_cmd = ' '.join(sort_cmd)

        with open(os.path.join(func_dir, 'sort.log'), 'w') as sort_log:
            subprocess.call(sort_cmd, shell=True, stderr=sort_log, stdout=sort_log)

        parsing_dir = os.path.join(m_config['ANNOTATION_DB_DIR'], 'parsing_db')

        eggnog_db = glob.glob(os.path.join(parsing_dir, '*eggnog*.db'))[0]
        refseq_db = glob.glob(os.path.join(parsing_dir, '*refseq*.db'))[0]
        refseq_tree = glob.glob(os.path.join(parsing_dir, '*.tree'))[0]

        out = os.path.join(blast_dir, os.path.splitext(db_name)[0] + '_parsed.tsv')
        if 'refseq' in db_name.lower():
            cmd = ['python', './scripts/blast2assignment.py', 'refseq', db_sorted, refseq_db, refseq_tree, out]

        elif 'eggnog' in db_name.lower():
            cmd = ['python', './scripts/blast2assignment.py', 'eggnog', db_sorted, eggnog_db, out]

        with open(os.path.join(blast_dir, 'blast_parsing_%s.log') % db_name, 'w') as flog:
            subprocess.call(cmd, stderr=flog, stdout=flog)


def merge_eggnog_refseq(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    func_dir = os.path.join(annotation_dir, 'functional')
    blast_dir = os.path.join(func_dir, 'blast_out')

    refseq_files = glob.glob(os.path.join(blast_dir, '*refseq*_parsed.tsv'))
    eggnog_files = glob.glob(os.path.join(blast_dir, '*eggnog*_parsed.tsv'))

    refseq_files.sort()
    eggnog_files.sort()

    for refseq, eggnog in zip(refseq_files, eggnog_files):
        out_name = read_pair_id + "_refseq_eggnog_merged.tsv"
        out = os.path.join(blast_dir, out_name)

        cmd = ['python', './scripts/blast2assignment.py', 'merge-tables', refseq, eggnog, out]
        with open(os.path.join(blast_dir, 'merge.log'), 'w') as f:
            subprocess.call(cmd, stderr=f, stdout=f)


def convert_gff_to_gtf(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    orf_dir = os.path.join(annotation_dir, 'ORFprediction')

    assembly_name = read_pair_id + '_final_contigs_1k.fasta'
    query = os.path.splitext(assembly_name)[0]

    gtf_name = query + '.gtf'
    gtf_path = os.path.join(orf_dir, gtf_name)

    gff_name = query + '.gff'
    gff_path = os.path.join(orf_dir, gff_name)

    inFile = gff_path
    outFile = gtf_path

    import re
    t = re.compile(r'ID=[0-9]+_([0-9]+);')
    with open(outFile, "w") as fh, open(inFile) as gff:
        print("##gff-version  3", file=fh)
        for line in gff:
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            orf = t.findall(toks[-1])[0]
            gene_id = toks[0] + "_" + orf
            toks[-1] = toks[-1] + "gene_id " + gene_id + ";"
            print(*toks, sep="\t", file=fh)


def verse(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    orf_dir = os.path.join(annotation_dir, 'ORFprediction')
    assembly_name = read_pair_id + '_final_contigs_1k.fasta'
    coverage_dir = get_coverage_dir(read_pair_id)

    q_dir = get_quant_dir(read_pair_id)
    if not os.path.exists(q_dir):
        os.mkdir(q_dir)

    query = os.path.splitext(assembly_name)[0]

    gtf_name = query + '.gtf'
    gtf_input = os.path.join(orf_dir, gtf_name)

    bam_input = glob.glob(os.path.join(coverage_dir, '*.bam'))

    numThreads = m_param['verse']['num_threads']
    min_read_overlap = m_param['verse']['min_read_overlap']

    with open(os.path.join(q_dir, 'verse.log'), 'w') as flog:
        for b in bam_input:
            out = os.path.join(q_dir, os.path.basename(b) + '.counts')
            cmd = ['verse', '--multithreadDecompress', '-T', numThreads, '--minReadOverlap', min_read_overlap,
                   '--singleEnd', '-t', 'CDS', '-z', '1', '-a', gtf_input, '-o', out, b]

            subprocess.call(cmd, stderr=flog, stdout=flog)


def funtaxa(read_pair_id):
    annotation_dir = get_annotation_dir(read_pair_id)
    func_dir = os.path.join(annotation_dir, 'functional')
    q_dir = get_quant_dir(read_pair_id)
    blast_dir = os.path.join(func_dir, 'blast_out')

    # define inputs
    count_tables = glob.glob(os.path.join(q_dir, '*.counts.CDS.txt'))
    annotation = os.path.join(blast_dir, read_pair_id + '_refseq_eggnog_merged.tsv')

    for count_table in count_tables:
        function_parser.get_counts(count_table, annotation)


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
    the_dir = os.path.join(the_dir, "Decon/")

    return the_dir


def get_flash_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Flash/")

    return the_dir


def get_trimmomatic_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Trimmomatic/")

    return the_dir


def get_fastqc_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "FastQC/")

    return the_dir


def get_interleave_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Interleave/")

    return the_dir


def get_assembly_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Assembly/")

    return the_dir


def get_coverage_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, 'Coverage/')

    return the_dir


def get_binning_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Binning/")

    return the_dir


def get_stats_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Stats/")

    return the_dir


def get_annotation_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Annotation/")

    return the_dir


def get_quant_dir(read_pair_id):
    the_dir = os.path.join(m_config['OUTPUT_DIR'], read_pair_id)
    the_dir = os.path.join(the_dir, "Quantification/")

    return the_dir


def get_database_dir():

    return os.path.abspath('./databases/')


def run_serial():
    performance_log = {}

    seperator = "-" * 100

    contaminant_dbs = m_param['decon']['contaminant_dbs']
    contaminant_dbs = contaminant_dbs.split(",")

    print("Building contaminant DBs...")
    for db_name in contaminant_dbs:
        db_name = db_name.strip()
        db_exists = contaminant_db_exists(db_name)
        if not db_exists:
            build_contaminant_db(db_name)

    print("Building annotation DBs...")
    build_annotation_dbs()

    for read_pair_id in m_read_pairs.keys():
        performance_log[read_pair_id] = {}

        print(seperator)
        print("Working on read pair '" + read_pair_id + "'")
        print(seperator + "\n")

        # Merge reads
        print("Merging reads...")
        start_time = time.time()
        merge_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['merging'] = end_time - start_time

        # Decon reads
        print("Decontaminating reads...")
        start_time = time.time()
        decon(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['decon'] = end_time - start_time

        # Trim reads
        print("Trimming reads...")
        start_time = time.time()
        trim_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['trimming'] = end_time - start_time

        # QC reads
        print("QC reads...")
        start_time = time.time()
        qc_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['qc'] = end_time - start_time

        # Interleave reads
        print("Interleave reads...")
        start_time = time.time()
        interleave_reads(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['interleaving'] = end_time - start_time

        print("Merging trim outputs...")
        start_time = time.time()
        merge_trim_outputs(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['mergeTrim'] = end_time - start_time

        print("Running MegaHit...")
        start_time = time.time()
        megahit(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['MegaHit'] = end_time - start_time

        print("Subsampling...")
        start_time = time.time()
        subsample(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['subsampling'] = end_time - start_time

        print("Generating statistics...")
        start_time = time.time()
        get_stats(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['stats'] = end_time - start_time

        print("Building assembly dbs...")
        start_time = time.time()
        build_assembly_db(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['buildAssemblyDB'] = end_time - start_time

        print("Mapping assembly reads...")
        start_time = time.time()
        map_to_assembly(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['mapAssembly'] = end_time - start_time

        print("Running MaxBin...")
        start_time = time.time()
        maxbin(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['MaxBin'] = end_time - start_time

        print("Running Prodigal..")
        start_time = time.time()
        prodigal(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['Prodigal'] = end_time - start_time

        print('Running Diamond...')
        start_time = time.time()
        diamond(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['Diamond'] = end_time - start_time

        print('Running Sword...')
        start_time = time.time()
        sword(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['Sword'] = end_time - start_time

        print('Blast parsing...')
        start_time = time.time()
        blast_parsing(read_pair_id)
        merge_eggnog_refseq(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['Blast Parsing'] = end_time - start_time

        print('Converting GFF to GTF...')
        start_time = time.time()
        convert_gff_to_gtf(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['GFF conversion'] = end_time - start_time

        print('Running verse...')
        start_time = time.time()
        verse(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['Verse'] = end_time - start_time

        print('Running funtaxa...')
        start_time = time.time()
        funtaxa(read_pair_id)
        end_time = time.time()

        performance_log[read_pair_id]['Funtaxa'] = end_time - start_time


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

    run_serial()
    # run_parallel()


if __name__ == "__main__":
    main()
