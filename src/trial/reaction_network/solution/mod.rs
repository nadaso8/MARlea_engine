use std::{collections::HashMap, fmt::Display, ops::Add};
use rayon::iter::ParallelIterator;

use super::reaction::Reaction;

/// Tuple struct wrapper around name data for a species of DNA
#[derive(PartialEq, Eq, PartialOrd, Ord,  Hash, Clone, Debug)]
pub struct Name(pub String);

/// Tuple struct wrapper around count data for a species of DNA
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Debug)]
pub struct Count(pub u64);

/// A collection of species name and counts which may be mutated over time by having Reaction objects applied to it via the .apply method 
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Solution {
    pub species_counts: HashMap<Name, Count>,
}

impl Solution {
    /// Mutates the solution to reflect the effects of a reaction
    pub fn apply (&mut self, reaction: Reaction) {
        for reactant in reaction.get_reactants() {
            self.species_counts.entry(reactant.get_species_name().clone())
                .and_modify(|species_count|
                    species_count.0 -= reactant.get_coefficient().0
                );
        }

        for product in reaction.get_products() {
            self.species_counts.entry(product.get_species_name().clone())
                .and_modify(|species_count|
                    species_count.0 += product.get_coefficient().0
                );
        }
    }

    /// Validates a reaction by checking it's reactants against the current solution 
    /// 
    /// #### Note!
    /// This function may default to true if only some of the reactants are in the current solution
    pub fn validate (&self, reaction: &Reaction) -> bool {
        let mut reaction_possible = true;

        for reactant in reaction.get_reactants() {
            if let Some(current_count) = self.species_counts.get(reactant.get_species_name()) {
                if reactant.get_coefficient().0 > current_count.0 {
                    reaction_possible = false;
                    break;
                }
            }
        }
        
        return reaction_possible;
    }
}

impl std::hash::Hash for Solution {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        for entry in self.species_counts.iter() {
            entry.hash(state);
        }
    }
}

impl IntoIterator for Solution {
    type Item = (Name, Count);
    type IntoIter = std::collections::hash_map::IntoIter<Name, Count>;

    fn into_iter(self) -> Self::IntoIter {
        self.species_counts.into_iter()
    }
    
}

impl Add for Solution {
    type Output = Self;

    /// preforms summation between Solutions in paralell by getting an entry from 
    fn add(self, rhs: Self) -> Self::Output {
        use rayon::iter::IntoParallelIterator;
        return Self{
            species_counts:
            self.species_counts
            .into_par_iter()
            .map(
                |(name, count)| 
                // calculates an entry< Name, Count > pair with a new count == count A + count B assuming count B is 0 if an entry for name DNE 
                ( name, Count( count.0 + rhs.species_counts.get( &name ).get_or_insert( &Count(0) ).0 ) ) 
            )
            .collect()
        }
    }
}

impl Display for Solution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Format the ordered vector as a string
        let mut formatted_string = String::new();
        for (name, count)in self.clone().into_iter() {
            formatted_string.push_str(&format!("{},{},", name.0.to_string(), count.0.to_string()));
        }

        // Write the formatted string to the provided Formatter
        write!(f, "{}", formatted_string)
    }
}
