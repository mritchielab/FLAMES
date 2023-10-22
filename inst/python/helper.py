import numpy as np
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
import os
import sys
import cProfile
import pstats
import time

def reverse_complement(seq):
	'''
	Args: <str>
		queried seq
	Returns: <str>
		reverse_complement seq
	'''
	comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
					'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	letters = \
		[comp[base] if base in comp.keys() else base for base in seq]
	return ''.join(letters)[::-1]

def err_msg(msg):
	CRED = '\033[91m'
	CEND = '\033[0m'
	print(CRED + msg + CEND)	

def warning_msg(msg):
	CRED = '\033[93m'
	CEND = '\033[0m'
	print(CRED + msg + CEND)

def green_msg(msg):
    CRED = '\033[92m'
    CEND = '\033[0m'
    print(CRED + msg + CEND)

def sliding_window_sum(array, window) :
    cum = np.cumsum(array)  
    return cum[window:] - cum[:-window]

def sliding_window_mean(array, window) :
    cum = np.cumsum(array)  
    return (cum[window:] - cum[:-window]) / window


class param:
    def __init__(self, **kwargs):
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
    def add(self, attr_name, attr_val, overwrite = True):
        if attr_name in __dict__.keys() and not overwrite:
            pass
        else:    
            self.__dict__[attr_name] = attr_val
    def rm(self, attr_name):
        if attr_name in self.__dict__.keys():
            del self.__dict__[attr_name]
    def __str__(self):
        return str(self.__dict__)
    
    def check(self, attr_list, add_none = True, silent = True):
        """
        Check whether the attributes in a given list are present. 

        Parameters
        ----------
        attr_list : LIST
            list of strings of the attributes to check 
        add_none : BOOL
            whether or not to create the missed attributes with value None
        silent : BOOL
            always return True if silent = True
        Returns
        -------
        True/False if all attributes is present
        """
        try:
            assert isinstance(add_none, bool)
            assert isinstance(silent, bool)
        except (AttributeError, TypeError):
            raise AssertionError("add_none and silent must be bool variable.")
        
        check_res = True
        for attr in attr_list:
            if attr not in self.__dict__.keys():
                check_res = False
                self.__dict__[attr] = None
        return check_res if not silent else True
    
def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,
                           pbar=True, pbar_unit='Read',pbar_func = lambda *x: 1, 
                           pbar_format=None,
                           schduler = 'process', 
                           preserve_order = True,
                           *arg, **kwargs):
    """multiple processing or threading, 
    Note: this function is adapted from BLAZE (github.com/shimlab/BLAZE)
    Args:
        func: function to be run parallely
        iterator: input to the function in each process/thread
        n_process (int, optional): number of cores or threads. Defaults to mp.cpu_count()-1.
        pbar (bool, optional): Whether or not to output a progres bar. Defaults to True.
        pbar_unit (str, optional): Unit shown on the progress bar. Defaults to 'Read'.
        pbar_func (function, optional): Function to calculate the total length of the progress bar. Defaults to len.
        schduler (str, optional): 'process' or 'thread'. Defaults to 'process'.

    Yields:
        return type of the func: the yield the result in the order of submit
    """
    if isinstance(iterator, (list, tuple)):
        iter_len = len(iterator)
        iterator = iter(iterator)
    else:
        iter_len = None
    
    if pbar:
        _pbar = tqdm(unit=pbar_unit, desc='Processed', 
                     total=iter_len, bar_format=pbar_format)

    # single process

    if n_process == 1:
        # a fake future object to make the code compatible with multiprocessing
        class fake_future:
            def __init__(self, rst):
                self.rst = rst
            def result(self):
                return self.rst
        for i in iterator:
            rst = func(i, *arg, **kwargs)
            if pbar:
                _pbar.update(pbar_func(i))
            yield fake_future(rst)
        return

    elif schduler == 'process':
        n_process = min(n_process, mp.cpu_count())
        executor = concurrent.futures.ProcessPoolExecutor(n_process)
    elif schduler == 'thread':
        executor = concurrent.futures.ThreadPoolExecutor(n_process)
    else:
        green_msg('Error in multiprocessing_submit: schduler should be either process or thread', printit=True)
        sys.exit(1)
     
    # A dictionary which will contain the future object
    max_queue = n_process + 10
    futures = {}
    n_job_in_queue = 0

    if preserve_order:
        # make sure the result is yield in the order of submit.
        job_idx = 0
        job_to_yield = 0
        job_completed = {}
        while True:
            while n_job_in_queue < max_queue:
                i = next(iterator, None)
                if i is None:
                    break
                futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
                job_idx += 1
                n_job_in_queue += 1

            # will wait until as least one job finished
            # batch size as value, release the cpu as soon as one job finished
            job = next(as_completed(futures), None)
            # no more job  
            if job is not None:
                job_completed[futures[job][1]] = job, futures[job][0]
                del futures[job]
                #print(job_completed.keys())
                # check order
                if  job_to_yield in job_completed.keys():
                    n_job_in_queue -= 1
                    # update pregress bar based on batch size
                    if pbar:
                        _pbar.update(job_completed[job_to_yield][1])
                    yield job_completed[job_to_yield][0]
                    del job_completed[job_to_yield]
                    job_to_yield += 1
            # all jobs finished: yield complelted job in the submit order
            else:
                while len(job_completed):
                    if pbar:
                        _pbar.update(job_completed[job_to_yield][1])
                    yield job_completed[job_to_yield][0]
                    del job_completed[job_to_yield]
                    job_to_yield += 1
                break
    else:
        while True:
            while n_job_in_queue < max_queue:
                i = next(iterator, None)
                if i is None:
                    break
                futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),i)
                n_job_in_queue += 1
            
            # will wait until as least one job finished
            # batch size as value, release the cpu as soon as one job finished
            job = next(as_completed(futures), None)
            # no more job  
            if job is not None:
                n_job_in_queue -= 1
                if pbar:
                    _pbar.update(futures[job][0])
                yield job
                del futures[job]
            else:
                break
# get file with a certian extensions
def get_files(search_dir, extensions, recursive=True):
    files = []
    if recursive:
        for i in extensions:
            files.extend(Path(search_dir).rglob(i))
        return files
    else:
        for i in extensions:
            files.extend(Path(search_dir).glob(i))
        return files

# check file exist
def check_exist(file_list):
    exit_code = 0
    for fn in file_list:
        if not os.path.exists(fn):
            exit_code = 1
            err_msg(f'Error: can not find {fn}')
    if exit_code == 1:
        sys.exit()


class pysam_AlignmentSegment_pickible:
    """Copy all necessary information from pysam AlignemntSegment to make it pickible 
    for multiprocessing
    """
    def __init__(self, AlignmentSegment):
        self.is_unmapped = AlignmentSegment.is_unmapped
        self.reference_name = AlignmentSegment.reference_name
        self.reference_end = AlignmentSegment.reference_end
        self.reference_start = AlignmentSegment.reference_start
        self.query_name = AlignmentSegment.query_name
        try:
            self.AS = AlignmentSegment.get_tag("AS")
        except:
            self.AS = None
        self.mapping_quality = AlignmentSegment.mapping_quality
        self.read_length = AlignmentSegment.infer_read_length()
        self.query_alignment_length = AlignmentSegment.query_alignment_length
    def get_tag(self, tag):
        if tag == "AS":
            return self.AS
        else:
            raise AttributeError
            sys.exit(1)
    def infer_read_length(self):
        return self.read_length
        
def batch_iterator(iterator, batch_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(pysam_AlignmentSegment_pickible(entry))
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):
        yield batch


def read_chunk_generator(file_handle, chunck_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    lines = []
    i=0
    chunk_number = 0
    for line in file_handle:
        i += 1
        lines.append(line)
        if i == chunck_size:
            # print(f"chunk {chunk_number} is ready")
            yield lines
            chunk_number += 1
            lines = []
            i = 0
    if len(lines):
        yield lines

def profile(fn):
    def profiled_fn(*args, **kwargs):
        
        profiler = cProfile.Profile()
        profiler.enable()
        start_time = time.time()
        result = fn(*args, **kwargs)
        start_time = time.time()
        end_time = time.time()
        wall_clock_time = end_time - start_time
        profiler.create_stats()

        # Change the filename as needed
        profile_filename = f"{fn.__name__}_profile.cprof"
        profiler.dump_stats(profile_filename)

        # Calculate the cpu percentage
        cputime = pstats.Stats(profile_filename).total_tt
        print(f"CPU percentage: {cputime/wall_clock_time*100}")
        return result
    return profiled_fn