use std::clone;

use reaction_network::ReactionNetwork;
use self::reaction_network::solution::Solution;

pub mod reaction_network; 

#[derive(Debug, Clone, Copy)]
pub struct Step(pub usize);
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ID(pub usize);

pub enum TrialState {
    Processing(ID, Solution),
    Complete(ID, Solution),
} 
/// The runtime environment for a single trial. Once the object has been initialized 
/// the step method must be called on it repeadedly until it reaches a stable state. 
#[derive(Clone)]
pub struct Trial {
    reaction_network: ReactionNetwork,
    stability: Stability, 
    max_semi_stable_steps: usize,
    trial_id: ID,
}

impl <'trial_runtime> Trial {

    pub fn from(
        reaction_network: ReactionNetwork,
        max_semi_stable_steps: usize, 
        trial_id: usize,
    ) -> Self {
        
        Self {
            reaction_network,
            stability: Stability::Initial,
            max_semi_stable_steps,
            trial_id: ID(trial_id), 
        }
    }

    /// Simulates a single time step in the trial and updates stability based on the list of possible reactions. 
    pub fn step(&mut self) -> TrialState {
        match self.stability {
            Stability::Initial => {
                self.reaction_network.react();

                if self.reaction_network.get_possible_reactions().is_empty() {
                    self.stability = Stability::Stable;
                }

                else if self.reaction_network.get_possible_reactions().is_subset(&self.reaction_network.get_null_adjacent_reactions()) {
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
                return TrialState::Complete(self.trial_id, self.reaction_network.get_solution().clone());
            }
        }
        return TrialState ::Processing(self.trial_id, self.reaction_network.get_solution().clone())
    }

    pub fn get_solution(&self) -> &Solution {
        self.reaction_network.get_solution()
    }

    pub fn get_seed(&self) -> [u8; 32] {
        self.reaction_network.get_seed().unwrap()
    }
}

#[derive(Clone, Copy)]
enum Stability {
    Initial, 
    Unstable,
    SemiStable(usize),
    Stable,
}