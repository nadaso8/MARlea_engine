use reaction_network::ReactionNetwork;
use results::TrialResult;
use std::sync::mpsc::SyncSender;

pub mod reaction_network; 
pub mod results;

/// Object specifiying return granularity for a trial
///  - Minimal returns only end stable solutions
///  - Full returns data for each step  
#[derive(Clone)]
pub enum TrialReturn {
    Minimal(SyncSender<TrialResult>),
    Full(SyncSender<TrialResult>)
}

/// The runtime environment for a single trial. Once the object has been initialized 
/// the simulate method may be called on it in order to simulate single CRN instance.
pub struct Trial {
    reaction_network: ReactionNetwork,
    stability: Stability, 
    step_counter: usize,
    max_semi_stable_steps: usize,
    id: usize,
    trial_return: TrialReturn
}

impl <'trial_runtime> Trial {

    pub fn from(
        reaction_network: ReactionNetwork,
        max_semi_stable_steps: usize, 
        id: usize,
        trial_return: TrialReturn
    ) -> Self {
        
        Self {
            reaction_network,
            stability: Stability::Initial,
            step_counter: 0,
            max_semi_stable_steps,
            id, 
            trial_return
        }
    }

    /// Calls the step method and returns the requested data until the reaction network enters the stable state. 
    pub fn simulate (&mut self)  {
        match self.trial_return.clone() {
            
            // behavior if minimal return behavior selected
            TrialReturn::Minimal(return_sender) => {
                loop{
                    self.step();
                    if let Stability::Stable = self.stability {
                        return_sender.send(TrialResult::StableSolution(self.reaction_network.get_solution().clone(), self.id, self.step_counter))
                        .expect("Reciever thread for trial {} dropped\nShutting down...");
                        return;
                    }
                }
            }

            // behavior if full histogram return selected
            TrialReturn::Full(return_sender) => {
                loop{
                    self.step();
                    return_sender.send(TrialResult::IntermediateStep(self.reaction_network.get_solution().clone(), self.id, self.step_counter))
                        .expect("Reciever thread for trial {} dropped\nShutting down...");
                    if let Stability::Stable = self.stability {
                        return_sender.send(TrialResult::StableSolution(self.reaction_network.get_solution().clone(), self.id, self.step_counter))
                        .expect("Reciever thread for trial {} dropped\nShutting down...");
                        return;
                    }
                }
            }
        }
   
    }

    /// Simulates a single time step in the trial and updates stability based on the list of possible reactions. 
    fn step(&mut self) {
        match self.stability {
            Stability::Initial => {
                self.reaction_network.react();

                if self.reaction_network.get_possible_reactions().is_empty() {
                    self.stability = Stability::Stable;
                }

                else if self.reaction_network.get_possible_reactions().is_subset(self.reaction_network.get_null_adjacent_reactions()) {
                    self.stability = Stability::SemiStable(0);
                } 
                
                else {
                    self.stability = Stability::Unstable;
                }
            } 

            Stability::Unstable => {
                self.reaction_network.react();

                if self.reaction_network.get_possible_reactions().is_empty() {
                    self.stability = Stability::Stable;
                }

                else if self.reaction_network.get_possible_reactions().is_subset(self.reaction_network.get_null_adjacent_reactions()) {
                    self.stability = Stability::SemiStable(0);
                }

                else {
                    self.stability = Stability::Unstable;
                }
            }

            Stability::SemiStable(count) => {
                self.reaction_network.react();

                if self.reaction_network.get_possible_reactions().is_empty() {
                    self.stability = Stability::Stable;


                } else if self.reaction_network.get_possible_reactions().is_subset(self.reaction_network.get_null_adjacent_reactions()) && count < self.max_semi_stable_steps {
                        self.reaction_network.react();
                        self.stability = Stability::SemiStable(count + 1);
                

                } else if self.reaction_network.get_possible_reactions().is_subset(self.reaction_network.get_null_adjacent_reactions()) && count >= self.max_semi_stable_steps {
                        self.reaction_network.react();
                        self.stability = Stability::Stable;
                

                } else {
                    self.stability = Stability::Unstable;
                }
            }

            Stability::Stable => {
                self.stability = Stability::Stable;
            }
        }
        self.step_counter += 1; 
    }
}

enum Stability {
    Initial, 
    Unstable,
    SemiStable(usize),
    Stable,
}