#!/user/bin/python2.7


# --- Imports ---
import sys, glob, shutil
import os.path

# --- Parameters ---
FLEXBAR_BARCODE_PARAM = ' -b /home/vincentm/aurelie/barcodes.fasta --threads 8 --barcode-unassigned --min-read-length 200 -be LEFT_TAIL -qt 22 --barcode-threshold 0 -t '
FLEXBAR_PRIMERS_PARAM = ' -b /home/vincentm/aurelie/primers.fasta --threads 8 --barcode-unassigned --min-read-length 200 -be LEFT -bk -qt 22 --barcode-threshold 0 --fasta-output -t '

# --- Constants ---
RFILE = 'Input_Fastq'
DEMULTIPLEX_BARCODE_DIR_PATH = 'Demultiplex_Barcode_Dir_Path'
DEMULTIPLEX_PRIMERS_DIR_PATH = 'Demultiplex_Primers_Dir_Path'
PRODUCTS = 'Products'
FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH = 'Flexbar_Out_Barcode_Fastq_File_Path'
FLEXBAR_OUT_BARCODE_FASTQ_FILE_NAME = 'FlexbarOut_barcode_Barcode_'
FLEXBAR_OUT_PRIMERS_RAW_FILE_NAME = ['FlexbarOut_Primers_barcode_ITS1_F.fasta', 'FlexbarOut_Primers_barcode_ITS2_F.fasta']
FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH = 'Flexbar_Out_Primers_Fasta_File_Path'
FLEXBAR_OUT_PRIMERS_FASTA_FILE_NAME = 'FlexbarOut_Primers_'



# --- Variables ---
data = dict()


for r_file in glob.glob('*.fastq'):
  ngs = int(r_file[-7:-6])
  data[ngs] = dict()
  # --- Create NGS folders ---
  data[ngs][DEMULTIPLEX_BARCODE_DIR_PATH] = os.path.join(os.getcwd(), 'NGS' + str(ngs) + '_Demultiplex_Barcode')
  if not os.path.exists(data[ngs][DEMULTIPLEX_BARCODE_DIR_PATH]):
    os.makedirs(data[ngs][DEMULTIPLEX_BARCODE_DIR_PATH])
  # --- RUN FLEXBAR-BARCODE ---
  # os.system('flexbar -r ' + \
            # os.path.join(os.getcwd() , r_file) + \
            # FLEXBAR_BARCODE_PARAM + \
            # os.path.join(data[ngs][DEMULTIPLEX_BARCODE_DIR_PATH], 'FlexbarOut'))

  data[ngs][PRODUCTS] = dict()
  for product in xrange(1, len(os.listdir(data[ngs][DEMULTIPLEX_BARCODE_DIR_PATH]))):
    data[ngs][PRODUCTS][product] = dict()
    data[ngs][PRODUCTS][product][FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH] = os.path.join(data[ngs][DEMULTIPLEX_BARCODE_DIR_PATH], \
                                                                                     FLEXBAR_OUT_BARCODE_FASTQ_FILE_NAME + str(product) + '.fastq')
    data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH] = os.path.join(os.getcwd(), 'NGS' + str(ngs) + '_Demultiplex_Primers')
    if not os.path.exists(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH]):
      os.makedirs(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH])

    # --- RUN FLEXBAR-PRIMERS ---
    
    # !-!-!  add if exist on fastq, else 'continue' !-!-!
    os.system('flexbar -r ' + \
              data[ngs][PRODUCTS][product][FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH] + \
              FLEXBAR_PRIMERS_PARAM + \
              os.path.join(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], 'FlexbarOut_Primers'))
    try:
      os.rename(os.path.join(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], 'FlexbarOut_Primers_barcode_unassigned.fasta'), \
                os.path.join(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], 'FlexbarOut_Primers_Unassigned_' + str(product) + '.fasta'))
    except OSError:
      pass
    
    for its in xrange(1, 3):
      data[ngs][PRODUCTS][product][its] = dict()
      data[ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH] = os.path.join(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], \
                                                                                            FLEXBAR_OUT_PRIMERS_FASTA_FILE_NAME + 'ITS' + str(its) + '_' + str(product) + '.fasta')
      try:
        os.rename(os.path.join(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], FLEXBAR_OUT_PRIMERS_RAW_FILE_NAME[its-1]), \
                  data[ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH])
      except OSError:
        print('###############################################')
        print('####### O S   E R R O R ######## F U C K ######')
        print('###############################################')
        print(os.path.join(data[ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], FLEXBAR_OUT_PRIMERS_RAW_FILE_NAME[its-1]))
        print('')
        print(data[ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH])
        print('###############################################')
        print('')
        print('')
        data[ngs][PRODUCTS][product].pop(its, None)
      
      # --- RELABEL ACCESS NUMBER ---
      try:
        relabeled_filename = os.path.join(os.path.split(data[ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH])[0], \
                                          'relabeled_NGS' + str(ngs) + '_ITS' + str(its) + '_' + os.path.basename(data[ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH]))
        os.system('vsearch ' + \
                  '--sortbylength ' + \
                  data[ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH] + \
                  ' --relabel NGS' + str(ngs) + '_ITS' + str(its) + '_ ' + \
                  '--maxseqlength 1000000 ' + \
                  '--output ' + relabeled_filename + \
                  ' --relabel_keep')
        # --- LINEARIZE FASTA FILE --- 
        linearised_filename = os.path.join(os.path.split(relabeled_filename)[0], \
                                           'linearised_' + os.path.basename(relabeled_filename))
        os.system('awk \'NR==1 {print ; next} {printf /^>/ ? "\\n"$0"\\n" : $1} END {printf "\\n"}\' ' + relabeled_filename + ' > ' + linearised_filename)
        # print('awk \'NR==1 {print ; next} {printf /^>/ ? "\\n"$0"\\n" : $1} END {printf "\\n"}\' ' + relabeled_filename + ' > ' + linearised_filename)
      
      except KeyError: 
        print('WARNING: No primers output for product ' + str(product) + ', NGS' + str(ngs) + ' ITS' + str(its) + '! ')

