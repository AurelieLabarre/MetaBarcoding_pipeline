#!/user/bin/python2.7


# --- Imports ---
import sys, glob, shutil
import os.path

# --- Parameters ---
# Leave space characters at beginning and end of parameter string!
FLEXBAR_BARCODE_PARAM = ' -b /home/vincentm/aurelie/barcodes.fasta --threads 8 --barcode-unassigned --min-read-length 200 -be LEFT_TAIL -qt 22 --barcode-threshold 0 -t '
FLEXBAR_PRIMERS_PARAM = ' -b /home/vincentm/aurelie/primers.fasta --threads 8 --barcode-unassigned --min-read-length 200 -be LEFT -bk -qt 22 --barcode-threshold 0 --fasta-output -t '
VSEARCH_DEREPLICATE_PARAM = ' --threads 8 --dbmask none --qmask none --rowlen 0 --notrunclabels --userfields query+id1+target --maxaccepts 0 --maxrejects 32 --top_hits_only --output_no_hits --db "/home/vincentm/aurelie/stampa/db_2_ITS_all.fasta" --id "0.5" --iddef 1 '
STAMPA_FOLDER = '/home/vincentm/aurelie/stampa'
PYTHON = 'python2.7 '

# --- Constants ---
RFILE = 'Input_Fastq'
DEMULTIPLEX_BARCODE_DIR_PATH = 'Demultiplex_Barcode_Dir_Path'
DEMULTIPLEX_PRIMERS_DIR_PATH = 'Demultiplex_Primers_Dir_Path'
FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH = 'Flexbar_Out_Barcode_Fastq_File_Path'
FLEXBAR_OUT_BARCODE_FASTQ_FILE_NAME = 'FlexbarOut_barcode_Barcode_'
FLEXBAR_OUT_PRIMERS_RAW_FILE_NAME = ['FlexbarOut_Primers_barcode_ITS1_F.fasta', 'FlexbarOut_Primers_barcode_ITS2_F.fasta']
FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH = 'Flexbar_Out_Primers_Fasta_File_Path'
FLEXBAR_OUT_PRIMERS_FASTA_FILE_NAME = 'FlexbarOut_Primers_'
CONCAT_DIR_PATH = 'Concat_Dir_Path'
SWARM_DIR_PATH = 'Swarm_Dir_Path'
SWARM_SORTED_DIR_PATH = 'Swarm_Sorted_Dir_Path'
PRODUCTS = 'Products'
LINEARISED_FULLPATH = 'Linearised_full_path'
NGS = 'NGS'


# --- Variables ---
data = dict()
data[NGS] = dict()
r_files = glob.glob('*.fastq')

for r_file in r_files:
  ngs = int(r_file[-7:-6])
  data[NGS][ngs] = dict()
  # --- Create NGS folders ---
  data[NGS][ngs][DEMULTIPLEX_BARCODE_DIR_PATH] = os.path.join(os.getcwd(), 'NGS' + str(ngs) + '_Demultiplex_Barcode')
  if not os.path.exists(data[NGS][ngs][DEMULTIPLEX_BARCODE_DIR_PATH]):
    os.makedirs(data[NGS][ngs][DEMULTIPLEX_BARCODE_DIR_PATH])
  # --- RUN FLEXBAR-BARCODE ---
  os.system('flexbar -r ' + \
            os.path.join(os.getcwd() , r_file) + \
            FLEXBAR_BARCODE_PARAM + \
            os.path.join(data[NGS][ngs][DEMULTIPLEX_BARCODE_DIR_PATH], 'FlexbarOut'))
  # --- Pause to modify barcode results ---
  wait = 0
  while(wait != ''):
    try:
      wait = raw_input('Press <Enter> to continue!')
    except:
      pass
  # ---
  data[NGS][ngs][PRODUCTS] = dict()
  for product in xrange(1, len(os.listdir(data[NGS][ngs][DEMULTIPLEX_BARCODE_DIR_PATH]))):
    data[NGS][ngs][PRODUCTS][product] = dict()
    data[NGS][ngs][PRODUCTS][product][FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH] = os.path.join(data[NGS][ngs][DEMULTIPLEX_BARCODE_DIR_PATH], \
                                                                                          FLEXBAR_OUT_BARCODE_FASTQ_FILE_NAME + str(product) + '.fastq')
    data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH] = os.path.join(os.getcwd(), 'NGS' + str(ngs) + '_Demultiplex_Primers')
    if not os.path.exists(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH]): 
      os.makedirs(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH])

    # --- RUN FLEXBAR-PRIMERS ---
    
    if os.path.exists(data[NGS][ngs][PRODUCTS][product][FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH]): 
      os.system('flexbar -r ' + \
                data[NGS][ngs][PRODUCTS][product][FLEXBAR_OUT_BARCODE_FASTQ_FILE_PATH] + \
                FLEXBAR_PRIMERS_PARAM + \
                os.path.join(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], 'FlexbarOut_Primers'))
    else: 
      continue
    try:
      os.rename(os.path.join(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], 'FlexbarOut_Primers_barcode_unassigned.fasta'), \
                os.path.join(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], 'FlexbarOut_Primers_Unassigned_' + str(product) + '.fasta'))
    except OSError:
      pass
    
    for its in xrange(1, 3):
      data[NGS][ngs][PRODUCTS][product][its] = dict()
      data[NGS][ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH] = os.path.join(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], \
                                                                                                 FLEXBAR_OUT_PRIMERS_FASTA_FILE_NAME + 'ITS' + str(its) + '_' + str(product) + '.fasta')
      try:
        os.rename(os.path.join(data[NGS][ngs][DEMULTIPLEX_PRIMERS_DIR_PATH], FLEXBAR_OUT_PRIMERS_RAW_FILE_NAME[its-1]), \
                  data[NGS][ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH])
      except OSError:
        data[NGS][ngs][PRODUCTS][product].pop(its, None)
      
      # --- RELABEL ACCESS NUMBER ---
      
      try:
        relabeled_filename = os.path.join(os.path.split(data[NGS][ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH])[0], \
                                          'relabeled_NGS' + str(ngs) + '_ITS' + str(its) + '_' + os.path.basename(data[NGS][ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH]))
        os.system('vsearch ' + \
                  '--sortbylength ' + \
                  data[NGS][ngs][PRODUCTS][product][its][FLEXBAR_OUT_PRIMERS_FASTA_FILE_PATH] + \
                  ' --relabel NGS' + str(ngs) + '_ITS' + str(its) + '_ ' + \
                  '--maxseqlength 1000000 ' + \
                  '--output ' + relabeled_filename + \
                  ' --relabel_keep')
        
        # --- LINEARIZE FASTA FILE --- 
        
        linearised_filename = os.path.join(os.path.split(relabeled_filename)[0], \
                                           'linearised_' + os.path.basename(relabeled_filename))
        os.system('awk \'NR==1 {print ; next} {printf /^>/ ? "\\n"$0"\\n" : $1} END {printf "\\n"}\' ' + relabeled_filename + ' > ' + linearised_filename)
        data[NGS][ngs][PRODUCTS][product][its][LINEARISED_FULLPATH] = linearised_filename

      except KeyError: 
        print('WARNING: No primers output for product ' + str(product) + ', NGS' + str(ngs) + ' ITS' + str(its) + '! ')
# --- Last NGS processed ---

# --- CONCATENATE NGS for each product - Keep ITS separated ---

max_products = 0
for ngs in data[NGS].keys(): 
  if len(data[NGS][ngs][PRODUCTS]) > max_products: 
    max_products = len(data[NGS][ngs][PRODUCTS])

data[CONCAT_DIR_PATH] = dict()
data[CONCAT_DIR_PATH][0] = os.path.join(os.getcwd(), 'Concatenated')
data[CONCAT_DIR_PATH][1] = os.path.join(data[CONCAT_DIR_PATH][0], 'ITS1')
data[CONCAT_DIR_PATH][2] = os.path.join(data[CONCAT_DIR_PATH][0], 'ITS2')
data[SWARM_DIR_PATH] = os.path.join(os.getcwd(), 'Swarm')
data[SWARM_SORTED_DIR_PATH] = os.path.join(os.getcwd(), 'Swarm_Sorted')
if not os.path.exists(data[CONCAT_DIR_PATH][0]):
  os.makedirs(data[CONCAT_DIR_PATH][0])
if not os.path.exists(data[CONCAT_DIR_PATH][1]):
  os.makedirs(data[CONCAT_DIR_PATH][1])
if not os.path.exists(data[CONCAT_DIR_PATH][2]):
  os.makedirs(data[CONCAT_DIR_PATH][2])
if not os.path.exists(data[SWARM_DIR_PATH]):
  os.makedirs(data[SWARM_DIR_PATH])
if not os.path.exists(data[SWARM_SORTED_DIR_PATH]):
  os.makedirs(data[SWARM_SORTED_DIR_PATH])

for product in xrange(1, max_products + 1): 
  for its in xrange(1, 3):
    concatenated_filename = os.path.join(data[CONCAT_DIR_PATH][its], \
                                         'concatenated_ITS' + str(its) + '_' + str(product) + '.fasta')
    fo = open(concatenated_filename, 'w')
    for ngs_cat in data[NGS].keys():
      try: 
        input = data[NGS][ngs_cat][PRODUCTS][product][its][LINEARISED_FULLPATH]
      except KeyError: 
        continue
      fi = open(input, 'r')
      lines = fi.readlines()
      fi.close()
      fo.writelines(lines)
    fo.close()

    # --- FORMING OTUs ---
    representative_filename = os.path.join(data[SWARM_DIR_PATH], \
                                           'representatives_' + os.path.basename(concatenated_filename))
    stat_filename = os.path.join(data[SWARM_DIR_PATH], \
                                 'stat_' + os.path.splitext(os.path.basename(concatenated_filename))[0] + '.txt')
    swarm_filename = os.path.join(data[SWARM_DIR_PATH], \
                                  'swarm_' + os.path.splitext(os.path.basename(concatenated_filename))[0] + '.txt')

    os.system('swarm -d 1 -f -w ' + representative_filename + ' -t 12 -s ' + stat_filename + ' < ' + concatenated_filename + ' > ' + swarm_filename)
        
        


    # --- Read swarm file into list ---
    input_file = swarm_filename
    f = open(input_file, 'r')
    # --- Delete empty lines ---
    lines = [line for line in f.readlines() if line.strip()]
    f.close()
    # --- Skip if empty file ---
    if lines == []: 
      continue

    singletons = []
    popped = 0
    for ii in xrange(len(lines)):
      # --- Create lists for each lines ---
      lines[ii-popped] = lines[ii-popped].split()
      if len(lines[ii-popped]) == 1:
        # --- Create list of singletons ---
        singletons.append(lines[ii-popped])
        # --- Delete singletons ---
        lines.pop(ii-popped)
        popped += 1
        # --- Progress ---
        if popped % 10000 == 0:
          sys.stdout.write('\rDeleting singletons: %d' % (popped))
          sys.stdout.flush()
    sys.stdout.write('\r')

    # --- Flatten singletons' list and add LF to each element ---
    singletons = [item + '\n' for sublist in singletons for item in sublist]
    # --- Write singletons' file ---
    output_singletons = os.path.join(data[SWARM_SORTED_DIR_PATH], 'singletons_ITS' + str(its) + '_' + str(product) + '.txt')
    f = open(output_singletons, 'w')
    f.writelines(singletons)
    f.close()
    sys.stdout.write('%d singletons deleted! Saved to %s\n' % (popped, output_singletons))

    ngs = []
    redondant_ngs = []
    popped = 0
    for ii in xrange(len(lines)):
      # --- Create list of lists of NGS ---
      ngs.append([lines[ii-popped][jj][3] for jj in xrange(len(lines[ii-popped]))])
      if len(set(ngs[ii])) < 2:
        # --- Create list of redondant NGS ---
        redondant_ngs.append(lines[ii-popped])
        # --- Delete redondant NGS ---
        lines.pop(ii-popped)
        popped += 1
        # --- Progress ---
        if popped % 1000 == 0:
          sys.stdout.write('\rDeleting redondants: %d' % (popped))
          sys.stdout.flush()
    sys.stdout.write('\r')

    for ii in xrange(len(redondant_ngs)):
      # --- Format redondant NGS data ---
      redondant_ngs[ii] = [redondant_ngs[ii][jj] + ' ' for jj in xrange(len(redondant_ngs[ii]))]
      redondant_ngs[ii][len(redondant_ngs[ii]) - 1] = redondant_ngs[ii][len(redondant_ngs[ii]) - 1].replace(' ', '\n')
    redondant_ngs = [item for sublist in redondant_ngs for item in sublist]
    # --- Write redondants' file ---
    output_redondants = os.path.join(data[SWARM_SORTED_DIR_PATH], 'redondants_ITS' + str(its) + '_' + str(product) + '.txt')
    f = open(output_redondants, 'w')
    f.writelines(redondant_ngs)
    f.close()
    sys.stdout.write('%d redondants deleted! Saved to %s\n' % (popped, output_redondants))

    # --- Find abondance of OTU ---
    abondance = [str(len(lines[ii])) + '\n' for ii in xrange(len(lines))]
    # --- Write to file ---
    abondance_file = os.path.join(data[SWARM_SORTED_DIR_PATH], 'abondance_ITS' + str(its) + '_' + str(product) + '.txt')
    f = open(abondance_file, 'w')
    f.writelines(abondance)
    f.close()
    sys.stdout.write('Abondance saved to %s\n' % (abondance_file))

    # --- Extract representative of OTU ---
    representative = [str(lines[ii][0]) + '\n' for ii in xrange(len(lines))]
    # --- Write to file ---
    representative_file = os.path.join(data[SWARM_SORTED_DIR_PATH], 'representatives_ITS' + str(its) + '_' + str(product) + '.txt')
    f = open(representative_file, 'w')
    f.writelines(representative)
    f.close()
    sys.stdout.write('Representatives saved to %s\n' % (representative_file))

    output_awesome = lines
    for ii in xrange(len(lines)):
      # --- Format output OTU without singletons and redondants ---
      output_awesome[ii] = [lines[ii][jj] + ' ' for jj in xrange(len(lines[ii]))]
      output_awesome[ii][len(lines[ii]) - 1] = output_awesome[ii][len(lines[ii]) - 1].replace(' ', '\n')
    output_awesome = [item for sublist in output_awesome for item in sublist]
    # --- Write output awesome file ---
    output_awesome_file = os.path.join(data[SWARM_SORTED_DIR_PATH], 'swarm_sorted_ITS' + str(its) + '_' + str(product) + '.txt')
    f = open(output_awesome_file, 'w')
    f.writelines(output_awesome)
    f.close()
    sys.stdout.write('Output awesome saved to %s\n' % (output_awesome_file))

    # --- Add abondance to representatives ---
    output = [lines[ii][0].strip() + '_' + str(len(lines[ii])) + '\n' for ii in range(len(lines))]
    output = [item for sublist in output for item in sublist]
    # --- Write output file ---
    output_file = os.path.join(data[SWARM_SORTED_DIR_PATH], 'representatives_abondance_ITS' + str(its) + '_' + str(product) + '.txt')
    f = open(output_file, 'w')
    f.writelines(output)
    f.close()
    sys.stdout.write('Representatives and their abondance saved to %s\n' % (output_file))

    # --- Open fasta file ---
    fasta_file = concatenated_filename
    f = open(fasta_file, 'r')
    fasta_lines = [line for line in f.readlines() if line.strip()]
    f.close()
    # --- Create dictionary with NGS as keys and lists of ID & sequence as values ---
    fasta = {}
    for ii in xrange(len(fasta_lines)): 
      fasta_lines[ii] = fasta_lines[ii].split()
      if ii % 2 == 0: 
        fasta[fasta_lines[ii][0][1:]] = [fasta_lines[ii][1], fasta_lines[ii+1].strip()]

    fasta_output = []
    for ii in xrange(len(lines)):
      fasta_output.append('>' + \
                          lines[ii][0].strip() + \
                          '_' + \
                          str(len(lines[ii])) + \
                          ' ')
      try: 
        fasta_output.append(fasta[lines[ii][0].strip()][0] + \
                            '\n' + \
                            fasta[lines[ii][0].strip()][1] + \
                            '\n')
      except KeyError:
        fasta_output.append("NO_ID_FOUND!\nNO_SEQUENCE_FOUND!\n")
        continue
    fasta_output = [item for sublist in fasta_output for item in sublist]
    output_fasta = os.path.join(data[SWARM_SORTED_DIR_PATH], 'final_to_stampa_ITS' + str(its) + '_' + str(product) + '.fasta')
    f = open(output_fasta, 'w')
    f.writelines(fasta_output)
    f.close()
    sys.stdout.write('Final fasta saved to %s\n' % (output_fasta))

    sys.stdout.flush()





    # --- DEREPLICATE (VSEARCH) ---
    assignments = os.path.join(os.path.split(output_fasta)[0], 'hits.' + os.path.split(output_fasta)[1])
    try:
      os.system('vsearch --usearch_global ' + \
                output_fasta + \
                VSEARCH_DEREPLICATE_PARAM + \
                '--userout ' + assignments)# + ' > "/dev/null" 2> "/dev/null"')
    except Exception as e: 
      print("Exception in vsearch bloc!")
      print(e.__doc__)
      print(e.message)
      pass
      # Since we are not splitting the file into chunks we don't need to do a merge
      # BUT the script does some parsing and sorting so we need to do that...
    try: 
      os.system(PYTHON + os.path.join(STAMPA_FOLDER, 'stampa_merge.py ') + data[SWARM_SORTED_DIR_PATH])
    except Exception as e: 
      print("Exception in python bloc!")
      print(e.__doc__)
      print(e.message)
      pass
      results = os.path.join(os.path.split(output_fasta)[0], 'results.' + os.path.split(output_fasta)[1])
    try:      
      os.system('sort --key=2,2nr --key=1,1d --output=' + results + ' ' + results)
    except Exception as e: 
      print("Exception in sort bloc!")
      print(e.__doc__)
      print(e.message)
      pass
