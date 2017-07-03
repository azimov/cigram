"""
get results for multiple samples of a given model with multiple processes
"""
import multiprocessing
import time


def sampler_function(M, m_args, measure_funcs, samples, results_dict, lock, seed):
    """
    Repeatedly run a model with a set of parameters and measure functions

    """
    for i in range(samples):
        t_start = time.time()
        m_args["seed"] = seed
        model_graph = M(**m_args)
        t_end = time.time()

        results = {}

        for measure, func in measure_funcs.items():
            results[measure] = func(model_graph)

        results["time"] = t_end - t_start

        if seed is not None:
            results["seed"] = seed
            seed +=1

        # store all the results at the same time point to keep list in the correct order
        lock.acquire()
        for measure, result in results.items():
            # multiprocess race condition needs lock
            # The following code looks odd but it breaks if you just do "results_dict[measure].append"
            # this is because of a weird bug in the multiprocessing manager shared dict class.
            # ommiting this will mean results_dict elements are lists of length 1.
            shared = results_dict[measure]
            shared.append(result)
            results_dict[measure]= shared
        lock.release()

        del model_graph


def multiprocess_samples(M, m_args, measure_funcs, samples=100, num_procs=multiprocessing.cpu_count(), seed=None):
    """
    run num_procs worth of samples
    M is the model function, m_args is the set of arugments to the model

    measure_funcs is the set of things we measure on the resulting graph object returned by the model

    if you want n_samples of graph e.g. do

    m_funcs = {
        "graphs":lambda G: G
    }

    this function would then return {"graphs":[<number of samples Graph objects>]}
    """
    if seed is None:
        seed = int(time())

    proc_sample_size = int(samples/num_procs)

    # more processes than samples
    if proc_sample_size == 0:
        proc_sample_size = 1

    manager = multiprocessing.Manager()
    results_dict = manager.dict()
    lock = manager.Lock()
    results_dict["time"] = []
    results_dict["seed"] = []

    for measure in measure_funcs:
        results_dict[measure] = []

    procs = []
    total_samples = 0
    for i in range(num_procs):
        # no need to spawn more processes than samples
        if total_samples >= samples:
            break

        n_samples = proc_sample_size
        total_samples += n_samples

        proc_seed = seed + (i * proc_sample_size)

        if i + 1 == num_procs:
            n_samples += samples - total_samples

        p = multiprocessing.Process(target=sampler_function, args=(M, m_args, measure_funcs, n_samples, results_dict, lock, proc_seed))
        procs.append(p)
        p.start()


    for p in procs:
        p.join()

    return dict(results_dict)
