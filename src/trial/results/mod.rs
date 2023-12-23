use super::reaction_network::solution::Solution;

type ID = usize;
type StepCounter = usize;

///The diferent types of data which may be sent back from a trial
#[derive(Eq, PartialEq, Clone)]
pub enum TrialResult {
    StableSolution(Solution, ID, StepCounter), 
    IntermediateStep(Solution, ID, StepCounter),
}