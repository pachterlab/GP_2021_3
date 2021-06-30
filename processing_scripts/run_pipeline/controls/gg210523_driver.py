import numpy as np
import multiprocessing
import sys
sys.path.append('../velocity-parameter-estimation/')
from seq_cme_inference import *
def param_parser(inputfile):
    with open(inputfile,'r') as f:
        L = f.readline()
        dataset_directory = read_split(f)
        result_directory = read_split(f)
        loom_filenames = read_split(f).split(',')
        transcriptome_filename = read_split(f)
        polyA_threshold = int(read_split(f))
        transcriptome_ind = int(read_split(f))
        filter_param = [float(k) for k in read_split(f).split(',')]
        all_prev_results = read_split(f).split(',')
        if len(all_prev_results[0])==0:
            all_prev_results=[]
        attribute_names = eval(read_split(f))
        gene_sel_seed = int(read_split(f))
        n_gen = int(read_split(f))
        loom_index = int(read_split(f))
        gene_result_list = read_split(f).split(',')
        if len(gene_result_list[0])==0:
            gene_result_list=[]
        phys_lb = np.array([float(k) for k in read_split(f).split(',')])
        phys_ub = np.array([float(k) for k in read_split(f).split(',')])
        search_restarts = int(read_split(f))
        init_pattern = read_split(f)
        use_lengths = bool(read_split(f))
        maxiter = int(read_split(f))
        n_pt1 = int(read_split(f))
        n_pt2 = int(read_split(f))
        samp_lb = np.array([float(k) for k in read_split(f).split(',')])
        samp_ub = np.array([float(k) for k in read_split(f).split(',')])
        ID_suffix = int(read_split(f))
        creator = read_split(f)
        NCOR = int(read_split(f))

    return dataset_directory, result_directory, loom_filenames,transcriptome_filename,\
    polyA_threshold, transcriptome_ind, filter_param, all_prev_results, attribute_names,\
    gene_sel_seed, n_gen, loom_index, gene_result_list, phys_lb, phys_ub, search_restarts,\
    init_pattern, use_lengths, maxiter, n_pt1,n_pt2,samp_lb, samp_ub, ID_suffix, creator, NCOR

def read_split(f):
	return f.readline().split(':')[1].strip()

def parout(PARINPUT):
	search_data,i=PARINPUT
	grid_search_driver(search_data,i)


def inference_workflow(input_param_file):
#	print(input_param_file)
	dataset_directory, result_directory, datasets, transcriptome_filename,\
    polyA_threshold, transcriptome_ind, filter_param, all_prev_results, attribute_names,\
    gene_sel_seed, n_gen, IND, gene_result_list, phys_lb, phys_ub, search_restarts,\
    init_pattern, use_lengths, maxiter,n_pt1,n_pt2, samp_lb, samp_ub, ID_suffix, creator, NCOR \
    = param_parser(input_param_file)
#	print(all_prev_results)
#	print(n_gen)
#	print(filter_param)
#	print(attribute_names)
	loom_filenames = [dataset_directory+k+'.loom' for k in datasets] #loom file
	print(loom_filenames)
#	print(IND)
#	print(gene_result_list)
	len_dict = get_transcriptome(transcriptome_filename)[transcriptome_ind]
	gene_set,trunc_gene_set = select_gene_set(loom_filenames,
	                      len_dict,viz=False,
	                          results_to_exclude=all_prev_results,seed=gene_sel_seed,n_gen=n_gen,
	                          filt_param=filter_param,attr_names_in=attribute_names)
	print('Gene set selected!')
	if len(gene_result_list)>0:
		result_data = import_datasets(gene_result_list)
		gene_set = result_data.gene_names
	search_data = get_gene_data(loom_filenames[IND],len_dict,gene_set,trunc_gene_set,viz=False,attr_names=attribute_names[IND])
	search_params = SearchParameters()
	search_params.define_search_parameters(search_restarts,phys_lb,phys_ub,maxiter,init_pattern,use_lengths)
	search_data.set_search_params(search_params)

	search_data.set_scan_grid(n_pt1,n_pt2,samp_lb,samp_ub)
	file_string = create_dir(search_data,result_directory,ID_suffix,meta=datasets[IND],creator=creator)

	search_data.get_pts()
	print('starting search...')
	t1 = time.time()
#	print(__name__)
#	if __name__ == 'inference_workflow':
#	if __name__ == '__main__':
	pool=multiprocessing.Pool(processes=NCOR)
	pool.map(parout,zip([search_data]*len(search_data.point_list),search_data.point_list))
	pool.close()
	pool.join()
	print('Parallelization done!')
	grid_search_driver(search_data)
	t2 = time.time()
	print('Runtime: {:.1f} min.'.format((t2-t1)/60))
	dump_results(file_string,include_nosamp=True)
	return

#print(__name__)
if __name__=="__main__":
#	print('hello')
#	print(__name__)
#	print(sys.argv)
	input_param_file=sys.argv[1]
#	print(input_param_file)
	inference_workflow(input_param_file)

