use std::marker::Send;
use std::sync::mpsc::Sender;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread::JoinHandle;

enum Message {
    NewJob(Job),
    Terminate,
}
use self::Message::{NewJob, Terminate};

pub struct ThreadPool {
    job_broadcast: Sender<Message>,
    join_handles: Vec<Option<JoinHandle<()>>>,
}

pub type Job = Box<FnBox + Send>;

pub trait FnBox {
    fn call_box(self: Box<Self>);
}

impl<F: FnOnce()> FnBox for F {
    fn call_box(self: Box<F>) {
        (*self)();
    }
}

impl ThreadPool {
    // create a thread pool with given thread count
    pub fn new(thread_count: usize) -> ThreadPool {
        assert!(thread_count > 0);

        let (job_broadcast, job_receive) = ::std::sync::mpsc::channel();

        let mut join_handles = vec![];

        let rx = Arc::new(Mutex::new(job_receive));

        for _i in 0..thread_count {
            let ri = Arc::clone(&rx);
            join_handles.push(Some(::std::thread::spawn(move || {
                loop {
                    let job = ri.lock().unwrap().recv().unwrap();

                    if let NewJob(job) = job {
                        // (*job)(); // this will compile one day, apparently
                        (job as Job).call_box();
                    } else {
                        break;
                    }
                }
            })));
        }

        ThreadPool {
            job_broadcast,
            join_handles,
        }
    }

    pub fn dispatch(&self, job: Job) {
        // broadcast a job to be picked up by the next thread
        self.job_broadcast
            .send(NewJob(job))
            .expect("Failed to broadcast job.");
    }
}

impl Drop for ThreadPool {
    fn drop(&mut self) {
        // signal termination to each thread by sending None as a job
        for _ in &self.join_handles {
            self.job_broadcast
                .send(Terminate)
                .expect("Failed to broadcast termination job.");
        }

        // come home little threadies
        for handle in &mut self.join_handles {
            if let Some(handle) = handle.take() {
                handle.join().expect("Failed to join thread");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn zero_thread_count_panics() {
        ThreadPool::new(0);
    }

    #[test]
    fn create_then_join() {
        ThreadPool::new(2);
    }

    #[test]
    fn run_trivial_job() {
        let job_ran_ref = Arc::new(Mutex::new(false));

        // dispatch a simple job, do a join, check that the job executed
        let tp = ThreadPool::new(2);

        {
            let job_ran_ref = Arc::clone(&job_ran_ref);

            tp.dispatch(Box::new(move || {
                // short sleep
                ::std::thread::sleep(::std::time::Duration::from_millis(20));
                *job_ran_ref.lock().unwrap() = true;
            }));
        }

        // job should not be finished yet due to sleep
        assert!(*job_ran_ref.lock().unwrap() == false);

        // wait for execution to finish and join thread
        drop(tp);

        // job flag should be set now
        assert!(*job_ran_ref.lock().unwrap());
    }
}
