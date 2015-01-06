import json

def convertJasparPFM2dict():
    """Parse Jaspar vertebrate pfm file and output to dict in JSON file

    location http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/
    """

    jaspar_file_name = "./pfm_vertebrates.txt"
    with open(jaspar_file_name, 'r') as jaspar_file_infile:
        jaspar_parse = jaspar_file_infile.readlines()    

    TFBS_matrix_dict = {}
    for i in range(0, len(jaspar_parse)):
        line = jaspar_parse[i]
        if ">" in line:
            TFBS_name = line[line.find(" ")+1:].replace('\n','')
            associated_matrix_strs = jaspar_parse[i+1:i+5]
            parsed_matrix = [[int(y) for y in x.split('\t')] for x in associated_matrix_strs]
            TFBS_matrix_dict[TFBS_name] = parsed_matrix

    updated_jaspar_file_name = "./TFBS_matrix.json"
    with open(updated_jaspar_file_name, 'w') as outfile:
        json.dump(TFBS_matrix_dict, outfile)
        

convertJasparPFM2dict()

