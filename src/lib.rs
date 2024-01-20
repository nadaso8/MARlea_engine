/// Author: Marceline Sorensen 
/// Email: nadaso8th@gmail.com 
/// Date: 12/21/2023
/// 
/// # Description
/// This Is the main simulation module for Marlea. It handles simulating to generate an estimated stable stat for CRNs via a monte carlo algorithm.
/// In order to use this module a ReactionNetwork object must first be constructed elsewhere before a runtime may be setup
/// 
/// # Known Issues 
/// When a marlea engine object is constructed solution data must contian all names in reaction reactant terms otherwise the reaction validation process will break.
/// 
/// # Hopeful Plans
/// I would like to switch the runtime of the backend over to rayon at some point 



pub mod trial;

use std::collections::{HashMap, HashSet};

use std::sync::mpsc::{
    sync_channel,
    SyncSender, Receiver
};
use std::{usize, thread};
use rayon::iter::ParallelIterator;
use rayon::{self, iter};
use trial::reaction_network::solution::{Name, Count, self};
use trial::{
    TrialState,
    reaction_network::{
        ReactionNetwork, 
        solution::Solution
    }
};

/// Marlea Error Types 
#[derive(Debug)]
pub enum MarleaError {
    Unknown(String)
}

/// The final average returned by a Marlea simulation
pub struct  MarleaResult(pub Vec<(String, f64)>);

/// The various behaviors marlea may use to return data as well any aditional data or objects they may need
#[derive(Clone)]
pub enum MarleaReturn {
    Full(SyncSender<MarleaResult>),
    Minimal(SyncSender<MarleaResult>),
    None()
}

/// This is a builder object containing defaults and methods for constructing a MarleaEngine Object. 
/// 
/// ## no response single threaded execution
/// 
/// ```
/// /* some code creating reactions and solution */
/// let reaction_network = ReactionNetwork::new(reactions, solution);
/// let marlea = marlea_engine::Builder::new(reaction_network).no_response()./*your options here*/.build();
/// 
/// /* you may drop the reciever since no_response option was selected causing all senders to also be dropped*/
/// drop(marlea.1);
/// 
/// /* simulate and return average of trials on current thread */
/// let result = marlea.0.run().unwrap().0;
/// 
/// ```
/// 
/// ## minimal - full response multithreaded execution
/// ```
/// /* some code creating reactions and solution */
/// let reaction_network = ReactionNetwork::new(reactions, solution);
/// let marlea = marlea_engine::Builder::new(reaction_network)./*your options here*/.build();
/// 
/// 
/// let frontend = threadpool::Threadpool::new(1);
/// frontend.execute( move|| {
///     /* my frontend code this should take the reciever offered up by the marlea_engine builder */
/// });
/// 
/// /* simulate and return average of trials on current thread */
/// let result = marlea.0.run().unwrap().0;
/// 
/// ```
pub struct Builder {
    // set externally
    num_trials: usize,
    max_runtime: Option<u64>,
    max_semi_stable_steps: usize,

    // constructed internally
    prime_network: ReactionNetwork,
    runtime_return: MarleaReturn,
    runtime_reciever: Receiver<MarleaResult>,
    trial_states: Vec<TrialState>,
}

impl Builder {
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
        let (runtime_sender, runtime_reciever) = sync_channel(128);
        let trial_states = Vec::new();

        Self { 
            num_trials: 100, 
            max_runtime: None, 
            max_semi_stable_steps: 99, 
            prime_network,
            runtime_return: MarleaReturn::Minimal(runtime_sender),
            runtime_reciever,
            trial_states,
        }
    }

    /// Sets the number of trials to simulate before taking an average and returning 
    pub fn trials(mut self, count: usize) -> Self {
        self.num_trials = count;
        self
    }
    
    /// Sets the maximum runtime to simulate for before forcing the simulation to quit
    pub fn runtime(mut self, time: u64) -> Self {
        self.max_runtime = Some(time);
        self
    }

    /// Sets the tolerance of stable step states to a manual value
    /// this should only be used in the case that a CRN is exiting prematurely and likely could be avoided with better network design
    pub fn tolerance(mut self, steps: usize) -> Self {
        self.max_semi_stable_steps = steps;
        self
    }

    /// Toggles ReturnVerbosity behavior between minimal and full
    pub fn verbose(mut self) -> Self {
        self.runtime_return = match self.runtime_return {
            MarleaReturn::Minimal(sender) => MarleaReturn::Full(sender),     
            MarleaReturn::Full(sender) => MarleaReturn::Minimal(sender),
            MarleaReturn::None() => MarleaReturn::None(),     
        };

        self
    }

    /// Sets Marlea to not offer any data back apart from the return from the .run() method 
    pub fn no_response(mut self) -> Self {
        match self.runtime_return.clone() {
            MarleaReturn::Full(sender) | MarleaReturn::Minimal(sender) => {
                drop(sender);
                self.runtime_return = MarleaReturn::None();
            },
            MarleaReturn::None() => ()
        }

        self
    }

    /// Consumes builder object and outputs a Marlea engine object 
    pub fn build(self) -> (MarleaEngine, Receiver<MarleaResult>) {

        let runtime = MarleaEngine {
            num_trials: self.num_trials,
            max_runtime:self.max_runtime,
            max_semi_stable_steps: self.max_semi_stable_steps,
            prime_network: self.prime_network,
            runtime_return: self.runtime_return,
            trial_states: self.trial_states
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
    prime_network: ReactionNetwork,
    runtime_return: MarleaReturn,
    trial_states: Vec<TrialState>,
}

impl MarleaEngine {
    /// Simulates the CRN and returns a sorted list of the average count of each species when a stable state is reached.  
    pub fn run(self) -> Result<MarleaResult, MarleaError>{

        // setup loop variables
        let mut trials_recieved = 0;
  
        // start runtime timer

        // setup trials 

        // step trials to completion with paralelliterator and monitor communication from frontend



        // attempt to return the final average 
        let result = MarleaResult(self.terminate());
        return Ok(result);
    }
    
    /// Averages the counts from a collection of solutions and returns them as a sorted list of name average value pairs 
    fn average_trials(sim_states: &Vec<Solution>) -> Vec<(String, f64)> {
        use rayon::iter::IntoParallelRefIterator;
        let summed_values = sim_states.par_iter().fold(||);
    }


    /// Cleanup runtime and print data
    fn terminate(self) -> Vec<(String, f64)> {
        
        let average_stable_solution = Self::average_trials(self.trial_states);

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

#[cfg(test)]
mod tests {
    use crate::trial::reaction_network::{reaction::{Reaction, term::Term}, solution::{self, Count, Name}};
    use super::*;

    #[test]
    fn sim_fibonacci_10 () {

        // manually described reaction network yes I know this is disgusting to look at but this should be able to run tests without influence from parser code
        let reactions = HashSet::from([
            Reaction::new( // csv ln 1
                vec![Term::new(Name("fibonacci.call".to_string()), Count(1)),],
                vec![Term::new(Name("setup.call".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 2
                vec![Term::new(Name("setup.done".to_string()), Count(1)),],
                vec![Term::new(Name("calculate.call".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 4
                vec![Term::new(Name("setup.call".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("next_value".to_string()), Count(1)),Term::new(Name("setup.call".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 5
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("next_value".to_string()), Count(2)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)),Term::new(Name("next_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 6
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("last_value".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 7
                vec![Term::new(Name("destruct".to_string()), Count(1)),Term::new(Name("current_value".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 8
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("setup.call".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 10
                vec![Term::new(Name("next_value.less_than.2.index.1".to_string()), Count(1)), Term::new(Name("setup.call.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("destruct.done.partial.0".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 11
                vec![Term::new(Name("next_value.less_than.2.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.less_than.2.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 12
                vec![Term::new(Name("next_value.less_than.2.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.less_than.2.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 13
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("next_value.less_than.2.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 14
                vec![Term::new(Name("next_value".to_string()), Count(2)), Term::new(Name("next_value.less_than.2.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("next_value".to_string()), Count(2)),],
                10000
            ),
            Reaction::new( // csv ln 15
                vec![Term::new(Name("next_value".to_string()), Count(2)), Term::new(Name("next_value.less_than.2.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("next_value".to_string()), Count(2)),],
                10000
            ),
            Reaction::new( // csv ln 16
                vec![Term::new(Name("setup.call.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("setup.call.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 17
                vec![Term::new(Name("setup.call.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("setup.call.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 18
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("setup.call.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 19
                vec![Term::new(Name("setup.call".to_string()), Count(1)), Term::new(Name("setup.call.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("setup.call".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 20
                vec![Term::new(Name("setup.call".to_string()), Count(1)),Term::new(Name("setup.call.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("setup.call".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 22    
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)), Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("destruct.done.partial.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 23
                vec![Term::new(Name("current_value.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 24
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 25
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("current_value.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 26
                vec![Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("current_value.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("current_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 27
                vec![Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("current_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 28
                vec![Term::new(Name("last_value.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 29
                vec![Term::new(Name("last_value.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 30
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("last_value.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 31
                vec![Term::new(Name("last_value".to_string()), Count(1)), Term::new(Name("last_value.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("last_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 32
                vec![Term::new(Name("last_value".to_string()), Count(1)), Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("last_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 34
                vec![Term::new(Name("destruct.done.partial.0".to_string()), Count(1)), Term::new(Name("destruct.done.partial.1".to_string()), Count(1)),],
                vec![Term::new(Name("destruct.done".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 35
                vec![Term::new(Name("destruct.done.partial.0".to_string()), Count(2)),],
                vec![Term::new(Name("destruct.done.partial.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 36
                vec![Term::new(Name("destruct.done.partial.1".to_string()), Count(2)),],
                vec![Term::new(Name("destruct.done.partial.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 37
                vec![Term::new(Name("destruct.done".to_string()), Count(2)),],
                vec![Term::new(Name("destruct.done".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 38
                vec![Term::new(Name("destruct.done".to_string()), Count(1)), Term::new(Name("destruct".to_string()), Count(1)),],
                vec![Term::new(Name("destruct.done".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 40
                vec![Term::new(Name("destruct.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("setup.done".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 41
                vec![Term::new(Name("destruct.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("destruct.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 42
                vec![Term::new(Name("destruct.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("destruct.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 43
                vec![Term::new(Name("destruct.done".to_string()), Count(1)),],
                vec![Term::new(Name("destruct.done".to_string()), Count(1)), Term::new(Name("destruct.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 44
                vec![Term::new(Name("destruct".to_string()), Count(1)),Term::new(Name("destruct.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 45
                vec![Term::new(Name("destruct".to_string()), Count(1)), Term::new(Name("destruct.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("destruct".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 46
                vec![Term::new(Name("setup.done".to_string()), Count(1)), Term::new(Name("destruct.done".to_string()), Count(1)),],
                vec![Term::new(Name("setup.done".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 48
                vec![Term::new(Name("calculate.call".to_string()), Count(2)),],
                vec![Term::new(Name("calculate.call".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 49
                vec![Term::new(Name("calculate.call".to_string()), Count(1)), Term::new(Name("calculate.done".to_string()), Count(1)),],
                vec![Term::new(Name("calculate.call".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 50
                vec![Term::new(Name("calculate.call".to_string()), Count(1)),],
                vec![Term::new(Name("index.check".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 51
                vec![Term::new(Name("index.check".to_string()), Count(1)),Term::new(Name("calculate.call".to_string()), Count(1)),],
                vec![Term::new(Name("index.check".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 53
                vec![Term::new(Name("index.check".to_string()), Count(2)),],
                vec![Term::new(Name("index.check".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 54
                vec![Term::new(Name("index.check".to_string()), Count(1)), Term::new(Name("index".to_string()), Count(1)),],
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 55
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)), Term::new(Name("index.check".to_string()), Count(1)),],
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 56
                vec![Term::new(Name("index.check".to_string()), Count(1)), Term::new(Name("index.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("calculate.return".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 57
                vec![Term::new(Name("index.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("index.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 58
                vec![Term::new(Name("index.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("index.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 59
                vec![Term::new(Name("index.check".to_string()), Count(1)),],
                vec![Term::new(Name("index.check".to_string()), Count(1)), Term::new(Name("index.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 60
                vec![Term::new(Name("index".to_string()), Count(1)), Term::new(Name("index.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("index".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 61
                vec![Term::new(Name("index".to_string()), Count(1)), Term::new(Name("index.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("index".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 62
                vec![Term::new(Name("calculate.return".to_string()), Count(1)), Term::new(Name("index.check".to_string()), Count(1)),],
                vec![Term::new(Name("calculate.return".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 64
                vec![Term::new(Name("current_value.convert".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 65
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)), Term::new(Name("current_value".to_string()), Count(1)),],
                vec![Term::new(Name("last_value".to_string()), Count(1)), Term::new(Name("current_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 66
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)), Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.convert".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 67
                vec![Term::new(Name("current_value.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 68
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 69
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)),],
                vec![Term::new(Name("current_value.convert".to_string()), Count(1)), Term::new(Name("current_value.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 70
                vec![Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("current_value.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("current_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 71
                vec![Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("current_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 73
                vec![Term::new(Name("next_value.convert".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 74
                vec![Term::new(Name("next_value.convert".to_string()), Count(1)), Term::new(Name("next_value".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.swap".to_string()), Count(1)), Term::new(Name("next_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 75
                vec![Term::new(Name("next_value.convert".to_string()), Count(1)), Term::new(Name("next_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.split".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 76
                vec![Term::new(Name("next_value.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 77
                vec![Term::new(Name("next_value.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 78
                vec![Term::new(Name("next_value.convert".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.not.index.0".to_string()), Count(1)), Term::new(Name("next_value.convert".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 79
                vec![Term::new(Name("next_value".to_string()), Count(1)), Term::new(Name("next_value.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("next_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 80
                vec![Term::new(Name("next_value".to_string()), Count(1)), Term::new(Name("next_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("next_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 81
                vec![Term::new(Name("next_value.split".to_string()), Count(1)), Term::new(Name("next_value.convert".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.split".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 83
                vec![Term::new(Name("next_value.split".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.split".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 84
                vec![Term::new(Name("next_value.split".to_string()), Count(1)), Term::new(Name("next_value.swap".to_string()), Count(1)),],
                vec![Term::new(Name("next_value".to_string()), Count(1)), Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("next_value.split".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 85
                vec![Term::new(Name("next_value.split".to_string()), Count(1)), Term::new(Name("next_value.swap.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 86
                vec![Term::new(Name("next_value.swap.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.swap.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 87
                vec![Term::new(Name("next_value.swap.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("next_value.swap.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 88
                vec![Term::new(Name("next_value.split".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.swap.not.index.0".to_string()), Count(1)), Term::new(Name("next_value.split".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 89
                vec![Term::new(Name("next_value.swap".to_string()), Count(1)), Term::new(Name("next_value.swap.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.swap".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 90
                vec![Term::new(Name("next_value.swap".to_string()), Count(1)), Term::new(Name("next_value.swap.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("next_value.swap".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 91
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)), Term::new(Name("next_value.split".to_string()), Count(1)),],
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 93
                vec![Term::new(Name("last_value.convert".to_string()), Count(2)),],
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 94
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)), Term::new(Name("last_value".to_string()), Count(1)),],
                vec![Term::new(Name("next_value".to_string()), Count(1)), Term::new(Name("last_value.convert".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 95
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)), Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("index.check".to_string()), Count(1)), Term::new(Name("last_value.convert".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 96
                vec![Term::new(Name("last_value.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 97
                vec![Term::new(Name("last_value.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 98
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)),],
                vec![Term::new(Name("last_value.convert".to_string()), Count(1)), Term::new(Name("last_value.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 99
                vec![Term::new(Name("last_value".to_string()), Count(1)), Term::new(Name("last_value.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("last_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 100
                vec![Term::new(Name("last_value".to_string()), Count(1)), Term::new(Name("last_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("last_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 101
                vec![Term::new(Name("index.check".to_string()), Count(1)), Term::new(Name("last_value.convert".to_string()), Count(1)),],
                vec![Term::new(Name("index.check".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 103
                vec![Term::new(Name("calculate.return".to_string()), Count(2)),],
                vec![Term::new(Name("calculate.return".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 104
                vec![Term::new(Name("calculate.return".to_string()), Count(1)), Term::new(Name("current_value".to_string()), Count(1)),],
                vec![Term::new(Name("return".to_string()), Count(1)), Term::new(Name("calculate.return".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 105
                vec![Term::new(Name("calculate.return".to_string()), Count(1)), Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("calculate.done".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 106
                vec![Term::new(Name("current_value.not.index.0".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 107
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(2)),],
                vec![Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 108
                vec![Term::new(Name("calculate.return".to_string()), Count(1)),],
                vec![Term::new(Name("calculate.return".to_string()), Count(1)), Term::new(Name("current_value.not.index.0".to_string()), Count(1)),],
                1
            ),
            Reaction::new( // csv ln 109
                vec![Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("current_value.not.index.0".to_string()), Count(1)),],
                vec![Term::new(Name("current_value".to_string()), Count(1)),],
                10000
            ),
            Reaction::new( // csv ln 110
                vec![Term::new(Name("current_value".to_string()), Count(1)), Term::new(Name("current_value.not.index.1".to_string()), Count(1)),],
                vec![Term::new(Name("current_value".to_string()), Count(1)),],
                10000
            ),
        ]);

        // manually described initial solution 
        let solution = Solution{species_counts: HashMap::from([
            (solution::Name("fibonacci.call".to_string()), Count(1)),
            (solution::Name("index".to_string()), Count(10)),
            (solution::Name("setup.call".to_string()), Count(0)),
            (solution::Name("setup.done".to_string()), Count(0)),
            (solution::Name("calculate.call".to_string()), Count(0)),
            (solution::Name("destruct".to_string()), Count(0)),
            (solution::Name("next_value".to_string()), Count(0)),
            (solution::Name("last_value".to_string()), Count(0)),
            (solution::Name("current_value".to_string()), Count(0)),
            (solution::Name("next_value.less_than.2.index.1".to_string()), Count(0)),
            (solution::Name("next_value.less_than.2.index.0".to_string()), Count(0)),
            (solution::Name("setup.call.not.index.1".to_string()), Count(0)),
            (solution::Name("setup.call.not.index.0".to_string()), Count(0)),
            (solution::Name("destruct.done.partial.1".to_string()), Count(0)),
            (solution::Name("destruct.done.partial.0".to_string()), Count(0)),
            (solution::Name("current_value.not.index.1".to_string()), Count(0)),
            (solution::Name("current_value.not.index.0".to_string()), Count(0)),
            (solution::Name("last_value.not.index.1".to_string()), Count(0)),
            (solution::Name("last_value.not.index.0".to_string()), Count(0)),
            (solution::Name("destruct.not.index.1".to_string()), Count(0)),
            (solution::Name("destruct.not.index.0".to_string()), Count(0)),
            (solution::Name("destruct.done".to_string()), Count(0)),
            (solution::Name("index.check".to_string()), Count(0)),
            (solution::Name("current_value.convert".to_string()), Count(0)),
            (solution::Name("index.not.index.1".to_string()), Count(0)),
            (solution::Name("index.not.index.0".to_string()), Count(0)),
            (solution::Name("calculate.return".to_string()), Count(0)),
            (solution::Name("calculate.done".to_string()), Count(0)),
            (solution::Name("next_value.convert".to_string()), Count(0)),
            (solution::Name("next_value.swap".to_string()), Count(0)),
            (solution::Name("next_value.split".to_string()), Count(0)),
            (solution::Name("next_value.not.index.1".to_string()), Count(0)),
            (solution::Name("next_value.not.index.0".to_string()), Count(0)),
            (solution::Name("next_value.swap.not.index.1".to_string()), Count(0)),
            (solution::Name("next_value.swap.not.index.0".to_string()), Count(0)),
            (solution::Name("last_value.convert".to_string()), Count(0)),
            (solution::Name("return".to_string()), Count(0)),
        ])};
        let (marlea, response_reciever) = Builder::new(ReactionNetwork::new(reactions, solution))
            .no_response()
            .trials(1000)
            .build();
        
        drop(response_reciever);
        
        let result = marlea.run().unwrap().0;
        
        for (name, count) in result {
            if name == "return".to_string() {
                assert!(count > 54.5&& count < 55.5);
            }
        }
    }
}