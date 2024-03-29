import multiprocessing
import threading
import time
from multiprocessing.managers import SyncManager
from multiprocessing.pool import ApplyResult, Pool
from multiprocessing.queues import Queue
from typing import Dict, Callable


def job_process(job_id, func: Callable, args_list: list, results_queue, stop_signal: multiprocessing.Value):
    results = list()
    for i, args in enumerate(args_list):
        # computation routine ...
        results.append(func(*args))

        # Update progress by sending it to the main process through a queue
        if (i + 1) % 100 == 0:
            results_queue.put((job_id, i + 1))

        if stop_signal.value:
            results_queue.put((job_id, 'stopped'))
            return list()

    results_queue.put((job_id, 'completed'))
    return results


class ProcessManager:
    def __init__(self, num_processes=5):
        self.manager: SyncManager = multiprocessing.Manager()
        self.results_queue: Queue = self.manager.Queue()
        self.pool: Pool = multiprocessing.Pool(num_processes)

        self.jobs: Dict[str, ApplyResult] = dict()
        self.jobs_stop_events: Dict[str, multiprocessing.Value] = dict()
        self.jobs_progress_callbacks: Dict[str, Callable] = dict()
        self.jobs_results_callbacks: Dict[str, Callable] = dict()
        self.jobs_args_list_length: Dict[str, int] = dict()

        self.monitor_progress_thread_stop = False
        self.monitor_progress_thread = threading.Thread(target=self._monitor_progress)
        self.monitor_progress_thread.start()

    def add_job(self, job_id, func, args_list, progress_callback, results_callback):
        if job_id in self.jobs.keys():
            return False

        stop_signal = self.manager.Value('i', 0)
        result = self.pool.apply_async(job_process, (job_id, func, args_list, self.results_queue, stop_signal))
        self.jobs[job_id] = result
        self.jobs_args_list_length[job_id] = len(args_list)
        self.jobs_stop_events[job_id] = stop_signal
        self.jobs_progress_callbacks[job_id] = progress_callback
        self.jobs_results_callbacks[job_id] = results_callback

        return True

    def _monitor_progress(self):
        while not self.monitor_progress_thread_stop:
            if self.results_queue.empty():
                time.sleep(0.1)
                continue

            job_id, progress = self.results_queue.get()

            if isinstance(progress, int):
                self.jobs_progress_callbacks[job_id](f"{job_id} progress: {progress}")
            if progress in ['completed', 'stopped']:
                self.jobs_results_callbacks[job_id](f'{job_id} results: {self.jobs[job_id].get()}')
                self._cleanup_job(job_id)

    def _cleanup_job(self, job_id: str):
        self.jobs.pop(job_id, None)
        self.jobs_args_list_length.pop(job_id, None)
        self.jobs_progress_callbacks.pop(job_id, None)
        self.jobs_results_callbacks.pop(job_id, None)
        self.jobs_stop_events.pop(job_id, None)

    def stop_job(self, job_id: str):
        self.jobs_stop_events[job_id].value = 1

    def shutdown(self):
        for job_id in list(self.jobs.keys()):
            self.stop_job(job_id)
        self.pool.close()
        self.pool.join()
        self.monitor_progress_thread_stop = True
        self.monitor_progress_thread.join()


def my_worker_func(v1):
    time.sleep(0.5)
    return v1


if __name__ == '__main__':
    pm = ProcessManager(num_processes=2)  # Limit to 2 concurrent processes for demonstration

    # Example usage
    args_list_ = [(i,) for i in range(5)]  # Dummy args for demonstration
    pm.add_job("Job 1", my_worker_func, args_list_, print, print)
    pm.add_job("Job 2", my_worker_func, args_list_, print, print)
    pm.add_job("Job 3", my_worker_func, args_list_, print, print)
    pm.add_job("Job 4", my_worker_func, args_list_, print, print)

    time.sleep(2)
    pm.stop_job('Job 2')
    pm.stop_job('Job 4')

    time.sleep(5)
    pm.add_job("Job 5", my_worker_func, args_list_, print, print)

    pm.shutdown()
