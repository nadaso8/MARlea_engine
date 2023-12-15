use super::Solution;

type ID = usize;
type StepCounter = usize;

#[derive(Eq, PartialEq, Clone)]
pub enum TrialResult {
    StableSolution(Solution, ID, StepCounter), 
    IntermediateStep(Solution, ID, StepCounter),
}