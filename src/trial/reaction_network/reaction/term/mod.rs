use super::super::solution::{Name, Count};

/// Contains The data for a single term within a reaction. 
#[derive(Debug, Hash, Eq, PartialEq, Clone, Ord, PartialOrd)]
pub struct  Term {
    species_name: Name,
    coefficient: Count,
}

impl Term {

    pub fn new(species_name: Name, coefficient: Count) -> Self {
        return Term{species_name , coefficient};
    }
    /// Returns a Count tuple struct reference
    pub fn get_coefficient (&self) -> &Count {
        return &self.coefficient;
    }

    /// Returns a Name tuple struct reference
    pub fn get_species_name(&self) -> &Name {
        return &self.species_name;
    }

}

