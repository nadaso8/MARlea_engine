/// Author: Marceline Sorensen 
/// Email: nadaso8th@gmail.com 
/// Date: 08/03/2023
/// 
/// # Description
/// This is the main simulation engine of MARlea it takes a set of reactions as well as a set of initial species valeus. 
/// It uses these to simulate the average stable case of a chemical reaction network stochastically. 
/// 
/// # Arguments (input_path, init_path)
/// - <input_path>
///     Specifies an input file location 
/// - <init_path> 
///     Specifies a file location to read initial species values from
///     Is of type Option which may be None
///     - If none will simply initialize all values to 0
/// - <out_path>
///     Specifies an output file location 
///     Is of type Option which may be None
/// - <num_trials>
///     Specifies the number of trials to be used in making a predicted average
///     Is of type Option which may be None
///     - if None will defualt to 100 trials to provide a simplistic estimation
/// - <max_runtime>
///     Specifies the maximum time the code my run for in seconds
///     Is of type Option which may be None 
///     - if None the simulation will run indefinitely
/// 
/// Accepted file types: 
///     - CSV
///     - XML UNIMPLEMENTED! 
///     - JSON UNIMPLEMENTED!


pub mod trial;

use std::collections::{HashMap, HashSet};

use std::sync::mpsc::{
    sync_channel,
    SyncSender, Receiver
};
use std::usize;
use threadpool::ThreadPool;
use trial::TrialReturn;
use trial::{
    results::TrialResult, 
    reaction_network::{
        ReactionNetwork, 
        reaction::
            term::solution::Solution
    }
};

/// Marlea Error Types 
pub enum MarleaError {
    Unknown(String)
}

type TrialID = usize;
type StepCounter = usize;
/// Marlea Response types  
#[derive(Clone)]
pub enum MarleaResponse {
    IntermediateStep(Solution, TrialID, StepCounter),
    StableSolution(Solution, TrialID, StepCounter),
    SimulationResult(Vec<(String, f64)>)
}

#[derive(Clone)]
pub enum MarleaReturn {
    Minimal(SyncSender<MarleaResponse>),
    Full(SyncSender<MarleaResponse>),
}

/// This is a builder object containing defaults and methods for constructing a MarleaEngine Object. 
/// 
/// # Usage 
/// ```MarleaEngineBuilder`
pub struct MarleaEngineBuilder {
    // set externally
    num_trials: usize,
    max_runtime: Option<u64>,
    max_semi_stable_steps: usize,

    // constructed internally
    computation_threads: ThreadPool,
    computation_threads_sender: SyncSender<TrialResult>,
    computation_threads_reciever: Receiver<TrialResult>,
    prime_network: ReactionNetwork,
    runtime_return: MarleaReturn,
    runtime_reciever: Receiver<MarleaResponse>
}

impl MarleaEngineBuilder {
    /// Builds a new MarleaEngine instance from given network with default values 
    /// 
    /// trials = 100
    /// runtime = unlimited 
    /// semi stable steps = 100
    /// return verbosity = minimal
    /// 
    pub fn new (
        prime_network: ReactionNetwork
    ) -> Self {
        let computation_threads = threadpool::Builder::new()
            .thread_name("compute_thread".to_string())
            .build();
        let (computation_threads_sender, computation_threads_reciever) = sync_channel(32);
        let (runtime_sender, runtime_reciever) = sync_channel(128);

        Self { 
            num_trials: 100, 
            max_runtime: None, 
            max_semi_stable_steps: 99, 
            computation_threads: computation_threads, 
            computation_threads_sender, 
            computation_threads_reciever, 
            prime_network,
            runtime_return: MarleaReturn::Minimal(runtime_sender),
            runtime_reciever,
        }
    }

    /// Sets the number of trials bto be executed to a manual value
    pub fn trials(mut self, count: usize) -> Self {
        self.num_trials = count;
        self
    }
    
    /// Sets the maximum runtime to a manual value
    pub fn runtime(mut self, time: u64) -> Self {
        self.max_runtime = Some(time);
        self
    }

    /// Sets the tolerance of stable step states to a manual value
    /// this should only be used in the case that a CRN is exiting prematurely 
    pub fn tolerance(mut self, steps: usize) -> Self {
        self.max_semi_stable_steps = steps;
        self
    }

    /// Toggles ReturnVerbosity behavior between minimal and full
    pub fn verbose(mut self) -> Self {
        self.runtime_return = match self.runtime_return {
            MarleaReturn::Minimal(sender) => MarleaReturn::Full(sender),     
            MarleaReturn::Full(sender) => MarleaReturn::Minimal(sender),     
        };

        self
    }

    /// Consumes builder object and outputs a Marlea engine object 
    pub fn build(self) -> (MarleaEngine, Receiver<MarleaResponse>) {

        let runtime = MarleaEngine {
            num_trials: self.num_trials,
            max_runtime:self.max_runtime,
            max_semi_stable_steps: self.max_semi_stable_steps,
            computation_threads: self.computation_threads,
            computation_threads_reciever: self.computation_threads_reciever,
            computation_threads_sender: self.computation_threads_sender,
            prime_network: self.prime_network,
            runtime_return: self.runtime_return
        };

        return(runtime, self.runtime_reciever)
    }
}

/// Main backend runtime object for Marlea
pub struct MarleaEngine {
    // set externally
    num_trials: usize,
    max_runtime: Option<u64>,
    max_semi_stable_steps: usize,

    // constructed internally
    computation_threads: ThreadPool,
    computation_threads_sender: SyncSender<TrialResult>,
    computation_threads_reciever: Receiver<TrialResult>,
    prime_network: ReactionNetwork,
    runtime_return: MarleaReturn,
}

impl MarleaEngine {
    /// Simulates the CRN and returns a sorted list of the average count of each species when a stable state is reached.  
    pub fn run(&self) -> Result<MarleaResponse, MarleaError>{

        // set containing all trial results
        let mut simulation_results = HashSet::new();

        // setup loop variables
        let mut trials_recieved = 0;
        let mut trials_created = 0;      
  
        // start runtime timer
        let (timer_sender, timer_reciever) = sync_channel(0);
        if let Some(time) = self.max_runtime {
            self.computation_threads.execute(move|| Self::engine_runtime_timer(time, timer_sender));
        }

        // setup trial return object
        let trial_return = match self.runtime_return {
            MarleaReturn::Minimal(_) => TrialReturn::Minimal(self.computation_threads_sender.clone()),
            MarleaReturn::Full(_) => TrialReturn::Full(self.computation_threads_sender.clone()),
        };
        while trials_created < self.num_trials {
            let mut current_trial = trial::Trial::from(
                self.prime_network.clone(), 
                self.max_semi_stable_steps,
                trials_created,
                trial_return.clone()
            );
            self.computation_threads.execute(move || current_trial.simulate());                    
            trials_created += 1;
        }

        // poll for trial results
        while trials_recieved < self.num_trials {
            if let Ok(result) = self.computation_threads_reciever.try_recv() {
                match &self.runtime_return {
                    MarleaReturn::Minimal(sender) => {
                        match result {
                            TrialResult::StableSolution(solution, id, step_counter) => {
                                trials_recieved += 1;
                                println!("Trial {} stable after {} steps", id, step_counter);
                                println!("Recieved {} trials", trials_recieved);
                                simulation_results.insert(solution.clone());
                                sender.send(MarleaResponse::StableSolution(solution, id, step_counter))
                                .expect("Fatal Error: frontend cannot be found by MarleaEngine. Please ensure the receiver is not dropped prematurely");
                            }
                            TrialResult::IntermediateStep(_,_,_) => (),
                        }
                    },
                    MarleaReturn::Full(sender) => {
                        match result {
                            TrialResult::StableSolution(solution, id, step_counter) => {
                                trials_recieved += 1;
                                println!("Trial {} stable after {} steps", id, step_counter);
                                println!("Recieved {} trials", trials_recieved);
                                simulation_results.insert(solution.clone());
                                sender.send(MarleaResponse::StableSolution(solution, id, step_counter))
                                .expect("Fatal Error: frontend cannot be found by MarleaEngine. Please ensure the receiver is not dropped prematurely");
                            }
                            TrialResult::IntermediateStep(solution, id, step_counter) => {
                                sender.send(MarleaResponse::IntermediateStep(solution, id, step_counter))
                                .expect("Fatal Error: frontend cannot be found by MarleaEngine. Please ensure the receiver is not dropped prematurely");
                            },
                        }
                    }, 
                }

            }
            
            if let Ok(_) = timer_reciever.try_recv() {
                println!("forced termination because max time was reached\n\nWARNING: returned results may not be accurate and should be used for debugging purposes only");
                break;
            }
        }

        // attempt to return the final average 
        let result = MarleaResponse::SimulationResult(self.terminate(simulation_results));
        match &self.runtime_return {
            MarleaReturn::Minimal(sender) | MarleaReturn::Full(sender) => {
                sender.send(result.clone())
                .expect("Fatal Error: frontend cannot be found by MarleaEngine. Please ensure the receiver is not dropped prematurely")
            },
            
        }
        return Ok(result);
    }
    
    fn average_trials(simulation_results: HashSet<Solution>) -> Vec<(String, f64)> {
        let mut summed_values = HashMap::<String, f64>::new();
        let num_trials = simulation_results.len() as f64;
    
        // Sum values of each species across all trials
        for result in &simulation_results {
            for (name, count) in &result.species_counts {
                summed_values.entry(name.clone())
                .and_modify(|summed_count| *summed_count += *count as f64)
                .or_insert(*count as f64);    
            }
        }
    
        // Calculate averages and sort alphabetically
        let mut averaged_values: Vec<(String, f64)> = summed_values.into_iter()
            .map(|(key, value)| (key, value / num_trials))
            .collect();
        averaged_values.sort_by_key(|(species, _)| species.to_owned());

        return averaged_values;
    }

    fn terminate(&self, simulation_results: HashSet<Solution>) -> Vec<(String, f64)> {
        
        let average_stable_solution = Self::average_trials(simulation_results);

        for (name, avg_count) in average_stable_solution.clone() {
            println!("{},{}", name , avg_count);
        }

        return average_stable_solution;
    }

    fn engine_runtime_timer(runtime: u64, tx: SyncSender<bool>) {
        let max_runtime = std::time::Duration::from_secs(runtime);
        std::thread::sleep(max_runtime);
        tx.send(true).unwrap();
        return;
    } 

}