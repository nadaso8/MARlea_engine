use super::super::solution::{
    Name,
    Count
};

/// Contains the data for a single term within a larger reaction.
/// Species is a reference to a named value in solution which will be added to or subtracted from. 
/// Coefficient is the value to add or subtract
#[derive(Debug, Hash, Eq, PartialEq, Clone)]
pub struct  Term {
    species_name: Name,
    coefficient: Count,
}

impl Term {

    /// creates a new term from a string slice
    /// returns none if the term sould be null
    /* DEPRICATED pub fn from(term: &str) -> Option<Self> {
        let mut species_name = None;
        let mut coefficient = None;
        let parts: Vec<&str> = term.split(" ").filter(|possible_part| !possible_part.is_empty()).collect();

        for part in parts {

            // if there is no coefficent try to parse as coefficent else parse as a name 
            let possible_coefficient = part.trim().parse::<u8>();
            let possible_name = part.trim().to_string();

            match possible_coefficient {
                Ok(value) => {
                    if let None = coefficient {coefficient = Some(value)}
                    else {panic!("more than one numeric value provided: it is unclear which is desired coefficient")}
                }
                Err(_) => {
                    if let None = species_name {species_name = Some(possible_name)}
                    else {
                        // non catastrophic error warn user with a print to console. 
                        print!("more than one possible name found in Term {}.\n{} was used as parsed name and coefficient was assumed to be 1.", term, species_name.clone().unwrap());
                    }
                }
            }

        }

        if let Some(parsed_name) = species_name {
            match coefficient {
                Some(parsed_value) => {
                    return Some(Term::new(parsed_name, parsed_value));
                },
                None => {
                    return Some(Term::new(parsed_name, 1 as u8));
                }
            }
        } else {return None}
    }*/

    pub fn new(species_name: Name, coefficient: Count) -> Self {
        return Term{species_name , coefficient};
    }
    /// returns the coefficient value of a Term
    pub fn get_coefficient (&self) -> &Count {
        return &self.coefficient;
    }

    /// Returns a reference to a Species enum
    /// Should always be a name
    pub fn get_species_name(&self) -> &Name {
        return &self.species_name;
    }

}

