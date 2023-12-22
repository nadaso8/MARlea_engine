use std::{collections::HashMap, fmt::Display};

use super::reaction::Reaction;

#[derive(PartialEq, Eq, PartialOrd, Ord,  Hash, Clone, Debug)]
pub struct Name(pub String);
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Debug)]
pub struct Count(pub u64);

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Solution {
    pub species_counts: HashMap<Name, Count>,
}

impl Solution {
    /// updates the counts in solution assuming the selected reaction was to happen
    pub fn apply (&mut self, reaction: Reaction) {
        for reactant in reaction.get_reactants() {
            self.species_counts.entry(reactant.get_species_name().clone())
                .and_modify(|species_count|
                    // if Counts aren't being updated it's likely that the internal u64 from the original count reference is being copied
                    // and not actually having it's value reassigned. 
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

    /// checks a reaction against the current state of solution and returns if it would be possible for the reaction to happen
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
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        // Make an ordered copy of self
        let mut self_as_vector: Vec<(Name, Count)> = Vec::new();
        for entry in self.species_counts {
            self_as_vector.push(entry);
        }
        self_as_vector.sort();

        return self_as_vector.into_iter();
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