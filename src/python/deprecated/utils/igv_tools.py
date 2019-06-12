from m6a.utils import m6a_utils
import csv

#PATH_MAPPING = '/Volumes/USBHUB/m6aAkker/Mapping/'
PATH_MAPPING = '/Volumes/m6aAkker/Mapping/'
#PATH_MAPPING = '/Users/UIBC/Documents/m6aExp/MappingBAM/'
#PATH_MAPPING = '/Volumes/SD/igv_m6a/BAM_file/'

def get_biocond_to_dataset_IPInput(exp_design_name):
    print(exp_design_name)
    file_name = m6a_utils.PATH + 'ExpDesign/' + exp_design_name + '_exp_design.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        biocond_to_dataset = dict()
        for row in exp_design:
            biocond = row.get('BioCond')
            sample_IP = row.get('IP_dataset')
            sample_Input = row.get('Input_dataset')
            sample = row.get('DataName')
            if biocond in biocond_to_dataset.keys():
                if not sample_ID in biocond_to_dataset[biocond]:
                    biocond_to_dataset[biocond].append(sample_ID)
            else:
                biocond_to_dataset[biocond] = [sample_ID]
        return biocond_to_dataset


def prepare_igv_all(exp_design_name):
    print(exp_design_name)
    biocond_to_dataset = dict()
    dataset_to_ip = dict()
    dataset_to_input = dict()
    file_name = m6a_utils.PATH + 'ExpDesign/' + exp_design_name + '_exp_design.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        for row in exp_design:
            biocond = row.get('BioCond')
            sample_IP = row.get('IP_dataset')
            sample_Input = row.get('Input_dataset')
            sample = row.get('DataName')
            if biocond in biocond_to_dataset.keys():
                if not sample in biocond_to_dataset[biocond]:
                    biocond_to_dataset[biocond].append(sample)
                    dataset_to_ip[sample].append(sample_IP)
                    dataset_to_input[sample].append(sample_Input)

            else:
                biocond_to_dataset[biocond] = [sample]
                dataset_to_ip[sample] = [sample_IP]
                dataset_to_input[sample] = [sample_Input]


    biocond_to_dataset = get_biocond_to_dataset_IPInput(exp_design_name)
    with open(m6a_utils.PATH_IGV + exp_design_name +'.xml', 'w') as igv_file:
        igv_file.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Global name=\"'+exp_design_name+'\" version=\"1\">\n')
        for biocond, datasets in biocond_to_dataset.items():

            igv_file.write('\t<Category name=\"'+biocond+'\">\n')
            igv_file.write('\t\t<Category name=\"Peaks\">\n')
            igv_file.write('\t\t\t<Resource name = \"' + biocond + ' Peaks\"\n\t\t\tpath = \"' + PATH_MAPPING +
                           'm6aExpDesign_' + exp_design_name + '_POI_'+biocond + '.bed\" />\n')
            igv_file.write('\t\t</Category>\n')
            for dataset in datasets:
                igv_file.write('\t\t<Category name=\"' + dataset + '\">\n')
                # igv_file.write('\t\t\t<Category name=\"POI score\">\n')
                # igv_file.write('\t\t\t\t<Resource name = \"'+dataset+' POI\"\n\t\t\t\tpath = \"'+PATH_MAPPING+
                #     dataset+'_POI.bedgraph\" />\n')
                # igv_file.write('\t\t\t</Category>\n')
                # igv_file.write('\t\t\t<Category name=\"POM score\">\n')
                # igv_file.write('\t\t\t\t<Resource name = \"' + dataset + ' POM IP\"\n\t\t\t\tpath = \"' + PATH_MAPPING +
                #                dataset + '_IP_POM.bedgraph\" />\n')
                # igv_file.write('\t\t\t\t<Resource name = \"' + dataset + ' POM Input\"\n\t\t\t\tpath = \"' + PATH_MAPPING +
                #                dataset + '_Input_POM.bedgraph\" />\n')
                # igv_file.write('\t\t\t</Category>\n')
                igv_file.write('\t\t\t<Category name=\"Norm. Coverage\">\n')
                igv_file.write('\t\t\t\t<Resource name = \"' + dataset + ' Cov IP\"\n\t\t\t\tpath = \"' + PATH_MAPPING +
                               dataset_to_ip.get(dataset) + '.bw\" />\n')
                igv_file.write('\t\t\t\t<Resource name = \"' + dataset + ' Cov Input\"\n\t\t\t\tpath = \"' + PATH_MAPPING +
                               dataset_to_input.get(dataset) + '.bw\" />\n')
                igv_file.write('\t\t\t</Category>\n')
                igv_file.write('\t\t\t<Category name=\"BAM file\">\n')
                igv_file.write('\t\t\t\t<Resource name = \"' + dataset + ' BAM IP\"\n\t\t\t\tpath = \"' + PATH_MAPPING +
                               dataset_to_ip.get(dataset) + '.bam\" />\n')
                igv_file.write('\t\t\t\t<Resource name = \"' + dataset + ' BAM Input\"\n\t\t\t\tpath = \"' + PATH_MAPPING +
                               dataset_to_input.get(dataset) + '.bam\" />\n')
                igv_file.write('\t\t\t</Category>\n')
                igv_file.write('\t\t</Category>\n')
            igv_file.write('\t</Category>\n')
        igv_file.write('</Global>')


def prepare_igv_specific(exp_design_name):
    print(exp_design_name)
    biocond_to_dataset = dict()
    dataset_to_ip = dict()
    dataset_to_input = dict()
    file_name = m6a_utils.PATH + 'ExpDesign/' + exp_design_name + '_exp_design.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        for row in exp_design:
            biocond = row.get('BioCond')
            sample_IP = row.get('IP_dataset')
            sample_Input = row.get('Input_dataset')
            sample = row.get('DataName')
            if biocond in biocond_to_dataset.keys():
                if not sample in biocond_to_dataset[biocond]:
                    biocond_to_dataset[biocond].append(sample)
                    dataset_to_ip[sample] = sample_IP
                    dataset_to_input[sample] = sample_Input

            else:
                biocond_to_dataset[biocond] = [sample]
                dataset_to_ip[sample] = sample_IP
                dataset_to_input[sample] = sample_Input

    # with open(m6a_utils.PATH_IGV + exp_design_name +'_POI.xml', 'w') as igv_file:
    #     igv_file.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Global name=\"'+exp_design_name+'_POI\" version=\"1\">\n')
    #     for biocond, datasets in biocond_to_dataset.items():
    #         igv_file.write('\t<Category name=\"'+biocond+'\">\n')
    #         for dataset in datasets:
    #             igv_file.write('\t\t<Category name=\"' + dataset + '\">\n')
    #             igv_file.write('\t\t\t<Resource name = \"'+dataset+' POI\"\n\t\t\tpath = \"'+PATH_MAPPING+
    #                 dataset+'_POI.bedgraph\" />\n')
    #             igv_file.write('\t\t</Category>\n')
    #         igv_file.write('\t</Category>\n')
    #     igv_file.write('</Global>')
    # with open(m6a_utils.PATH_IGV + exp_design_name +'_POM.xml', 'w') as igv_file:
    #     igv_file.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Global name=\"'+exp_design_name+'_POM\" version=\"1\">\n')
    #     for biocond, datasets in biocond_to_dataset.items():
    #         igv_file.write('\t<Category name=\"'+biocond+'\">\n')
    #         for dataset in datasets:
    #             igv_file.write('\t\t<Category name=\"' + dataset + '\">\n')
    #             igv_file.write('\t\t\t<Resource name = \"'+dataset+' POM IP\"\n\t\t\tpath = \"'+PATH_MAPPING+
    #                 dataset+'_IP_POM.bedgraph\" />\n')
    #             igv_file.write('\t\t\t<Resource name = \"' + dataset + ' POM Input\"\n\t\t\tpath = \"' + PATH_MAPPING +
    #                            dataset + '_Input_POM.bedgraph\" />\n')
    #             igv_file.write('\t\t</Category>\n')
    #         igv_file.write('\t</Category>\n')
    #     igv_file.write('</Global>')
    with open(m6a_utils.PATH_IGV + exp_design_name +'_Cov.xml', 'w') as igv_file:
        igv_file.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Global name=\"'+exp_design_name+'_Norm. Coverage\" version=\"1\">\n')
        for biocond, datasets in biocond_to_dataset.items():
            igv_file.write('\t<Category name=\"' + biocond + '\">\n')
            igv_file.write('\t\t<Category name=\"' + biocond + ' Mean\">\n')
            igv_file.write(
                    '\t\t\t<Resource name = \"' + biocond + ' Mean Norm Cov IP\"\n\t\t\tpath = \"' + PATH_MAPPING +
                    biocond + '_IP.bw\" />\n')
            igv_file.write(
                    '\t\t\t<Resource name = \"' + biocond + ' Mean Cov Input\"\n\t\t\tpath = \"' + PATH_MAPPING +
                    biocond + '_Input.bw\" />\n')
            igv_file.write('\t\t</Category>\n')
            igv_file.write('\t\t<Category name=\"' + biocond + ' SD\">\n')
            igv_file.write(
                '\t\t\t<Resource name = \"' + biocond + ' SD Norm Cov IP\"\n\t\t\tpath = \"' + PATH_MAPPING +
                biocond + '_IP_SD.bw\" />\n')
            igv_file.write(
                '\t\t\t<Resource name = \"' + biocond + ' SD Cov Input\"\n\t\t\tpath = \"' + PATH_MAPPING +
                biocond + '_Input_SD.bw\" />\n')
            igv_file.write('\t\t</Category>\n')
            igv_file.write('\t\t<Category name=\"' + biocond + ' Mean_Old\">\n')
            igv_file.write(
                '\t\t\t<Resource name = \"' + biocond + ' Mean Old Norm Cov IP\"\n\t\t\tpath = \"' + PATH_MAPPING +
                biocond + '_IP_Old.bw\" />\n')
            igv_file.write('\t\t</Category>\n')
            for dataset in datasets:
                igv_file.write('\t\t<Category name=\"' + dataset + '\">\n')
                igv_file.write('\t\t\t<Resource name = \"'+dataset+' Norm Cov IP\"\n\t\t\tpath = \"'+PATH_MAPPING+
                               dataset_to_ip.get(dataset)+'.bw\" />\n')
                igv_file.write('\t\t\t<Resource name = \"' + dataset + ' Norm Cov Input\"\n\t\t\tpath = \"' + PATH_MAPPING +
                               dataset_to_input.get(dataset) + '.bw\" />\n')
                igv_file.write('\t\t</Category>\n')
            igv_file.write('\t</Category>\n')
        igv_file.write('</Global>')
    with open(m6a_utils.PATH_IGV + exp_design_name +'_BAM.xml', 'w') as igv_file:
        igv_file.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Global name=\"'+exp_design_name+'_BAM\" version=\"1\">\n')
        for biocond, datasets in biocond_to_dataset.items():
            igv_file.write('\t<Category name=\"'+biocond+'\">\n')
            for dataset in datasets:
                igv_file.write('\t\t<Category name=\"' + dataset + '\">\n')
                igv_file.write('\t\t\t<Resource name = \"'+dataset+' BAM IP\"\n\t\t\tpath = \"'+PATH_MAPPING+
                               dataset_to_ip.get(dataset)+'.bam\" />\n')
                igv_file.write('\t\t\t<Resource name = \"' + dataset + ' BAM Input\"\n\t\t\tpath = \"' + PATH_MAPPING +
                               dataset_to_input.get(dataset) + '.bam\" />\n')
                igv_file.write('\t\t</Category>\n')
            igv_file.write('\t</Category>\n')
        igv_file.write('</Global>')
    # with open(m6a_utils.PATH_IGV + exp_design_name +'_Peaks.xml', 'w') as igv_file:
    #     igv_file.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Global name=\"'+exp_design_name+'_Peaks\" version=\"1\">\n')
    #     for biocond, datasets in biocond_to_dataset.items():
    #         igv_file.write('\t<Category name=\"'+biocond+'\">\n')
    #         igv_file.write('\t\t<Resource name = \"' + biocond + ' Peaks\"\n\t\tpath = \"' + PATH_MAPPING
    #                        + exp_design_name + '_POI_' + biocond + '.bed\" />\n')
    #         igv_file.write('\t</Category>\n')
    #     igv_file.write('</Global>')


#prepare_igv_all('Seq6')
#prepare_igv_specific('Seq6')
#prepare_igv_all('Seq789')
#prepare_igv_specific('Seq789')

#prepare_igv_specific('Cecum')
#prepare_igv_specific('Liver')
prepare_igv_specific('CecAm')
prepare_igv_specific('LivOld')
# prepare_igv_all('Liver_ZT3')
# prepare_igv_all('Liver_ZT11')
# prepare_igv_specific('Liver_ZT3')
# prepare_igv_specific('Liver_ZT11')