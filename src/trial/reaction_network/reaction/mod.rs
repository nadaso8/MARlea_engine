pub mod term;

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use term::Term;

/// Represents a single reaction, containing a list of reactants and products as well as a reaction rate
/// in order to scale the probability of that reaction occurring.
#[derive(Debug, Eq, PartialEq, Clone, Ord, PartialOrd)]
pub struct Reaction {
    reactants: Vec<Term>,
    products: Vec<Term>,
    reaction_rate: u64,
}

impl Reaction {

    pub fn new (reactants: Vec<Term>, products: Vec<Term>, reaction_rate: u64) -> Self {
        return Self { reactants: reactants, products: products, reaction_rate: reaction_rate};
    }
    
    /// Returns a reference to the list of reactants for a reaction
    pub fn get_reactants(&self) -> &Vec<Term> {
        return &self.reactants;
    }

    /// Returns a reference to the list of products for a reaction
    pub fn get_products(&self) -> &Vec<Term> {
        return &self.products;
    }

    /// Returns the reaction rate for a reaction
    pub fn get_reaction_rate (&self) -> u64 {
        return self.reaction_rate;
    }
}

impl Hash for Reaction {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut hasher = DefaultHasher::new();
        for term in &self.reactants {
            term.hash(&mut hasher);
        }
        for term in &self.products {
            term.hash(&mut hasher);
        }
        self.reaction_rate.hash(&mut hasher);
        hasher.finish().hash(state);
    }
}