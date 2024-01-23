use std::collections::{HashSet, BTreeSet};
use rand::{Rng, rngs::StdRng, SeedableRng};
use reaction::Reaction;
use solution::Solution;

pub mod solution;
pub mod reaction; 

/// Data structure representing a DNA chemical reaction network. 
/// - reactions 
///     - the set of reactions in the network
/// - possible reactions
///     - the set of reactions which may happen next
/// - null_adjacent reactions 
///     - the set of reactions which will always happen or contain the products of reactions which always happen.
///       if the possible reactions is a subset of this it indicates the reaction network may have stalled out to a stable state
#[derive(Clone)]
pub struct ReactionNetwork {
    reactions: HashSet<Reaction>, 
    possible_reactions: BTreeSet<Reaction>,
    null_adjacent_reactions: BTreeSet<Reaction>,
    solution: Solution,
    prng: StdRng,
    seed: [u8; 32],
}

impl ReactionNetwork {

    pub fn new(reactions: HashSet<Reaction>, solution: Solution) -> Self {
        // Initialize null_adjacent_reactions and possible_reactions as empty HashSet
        let null_adjacent_reactions = BTreeSet::new();
        let possible_reactions = BTreeSet::new();
        let seed: [u8; 32] = rand::random();
        let prng = rand::rngs::StdRng::from_seed(seed);

        // Make a new instance of Self with the provided arguments and initialized fields.
        let mut new_netowrk = 
        Self{
            reactions, 
            solution, 
            possible_reactions, 
            null_adjacent_reactions,
            prng,
            seed,
        };

        // Generate and cache null adjacent and possible reactions up front so they are always available
        new_netowrk.gen_null_adjacent_reactions();
        new_netowrk.find_possible_reactions();

        return new_netowrk;
    }


    /// sets prng to some prng based on seed
    pub fn with_seed(mut self, seed: [u8; 32]) -> Self {
        self.seed = seed;
        self.prng = StdRng::from_seed(seed);
        self
    }

    /// returns a reference to the current seed value
    pub fn get_seed(&self) -> &[u8; 32] {&self.seed}

    /// Returns a reference to the set of null adjacent reactions
    pub fn get_null_adjacent_reactions(&self) -> &BTreeSet<Reaction> {
        return &self.null_adjacent_reactions;
    }

    /// Determines and caches the set of null adjacent reactions
    fn gen_null_adjacent_reactions(&mut self) {

        self.null_adjacent_reactions.clear();

        for reaction in &self.reactions {
            // Check for reactions that only have products (null adjacent).
            if reaction.get_reactants().is_empty() {

                // Insert the reaction into the null_adjacent_reactions HashSet and access its corresponding product(s)
                if self.null_adjacent_reactions.insert(reaction.clone()) {
                    for product in reaction.get_products() {
                        let null_generated_species = product.get_species_name().0.clone();

                        // For each secondary reaction, check if its reactant species matches the current null generated species
                        for secondary_reaction in &self.reactions {
                            for secondary_reactant in secondary_reaction.get_reactants() {

                                if null_generated_species == secondary_reactant.get_species_name().0.clone() {
                                    // Insert the reaction into the null_adjacent_reactions HashSet.
                                    self.null_adjacent_reactions.insert(secondary_reaction.clone());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Returns the cached set of possible reactions
    pub fn get_possible_reactions (&self) -> &BTreeSet<Reaction>{
        &self.possible_reactions
    }

    /// Re-generates the list of possible reactions
    fn find_possible_reactions (&mut self) {
        self.possible_reactions.clear();

        // loop over all reactions and check if it's possible for them to occur based on current species concentration
        for reaction in &self.reactions {
            if self.solution.validate(reaction) {
                self.possible_reactions.insert(reaction.clone()); // add reaction to list of possible reactions
            }
        }
    }

    /// Takes the sum of the all the reaction rates currently possible_reactions 
    fn sum_possible_reaction_rates (&self) -> u128 {
        let mut sum: u128 = 0; 
        // loop over all possible reactions and sum their reaction rates
        for reaction in self.get_possible_reactions() {
            sum += reaction.get_reaction_rate() as u128;
        }
        return sum;
    }


    /// Randomly selects a reaction from possible_reactions to happen next weighted by the reaction rate of each reaction
    fn get_next_reaction (&mut self) -> Result<Reaction, String> {

        let max_index = self.sum_possible_reaction_rates();
        let mut index = self.prng.gen_range(0..max_index);
        let mut next_reaction: Result<Reaction, String> = Result::Err("failed to find next reaction".to_string());

        // iterate through all possible valid reactions and pick one based on its probability 
        for reaction in self.get_possible_reactions() {
            if reaction.get_reaction_rate() as u128 > index {
                next_reaction = Result::Ok(reaction.clone());
                break;
            } else {
                index -= reaction.get_reaction_rate() as u128;
            }
        }

        return next_reaction;
    }

    /// Selects a reaction to happen next, applies it to the solution, and re-generates the set of possible reactions so it remains up to date. 
    pub fn react(&mut self) {
        let next_reaction = self.get_next_reaction().unwrap();
        self.solution.apply(next_reaction);
        self.find_possible_reactions();// update list of possible reactions after solution has changed
    }

    // Returns a reference to the current solution
    pub fn get_solution(&self) -> &Solution {
        return &self.solution;
    }
}

