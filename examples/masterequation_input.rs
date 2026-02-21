use MarXus::masterequation::api::{run_multiwell_from_network, MultiWellFromNetworkInput};
use MarXus::masterequation::input_deck::parse_master_equation_input;
use MarXus::masterequation::reaction_network::ChemicalActivationDefinition;
use MarXus::masterequation::example::PlaceholderMicroData;

fn main() -> Result<(), String> {
    let input = include_str!("masterequation_input.txt");
    let network = parse_master_equation_input(input)?;

    let micro = PlaceholderMicroData;

    let activation = ChemicalActivationDefinition {
        activated_well_index: 0,
        recombination_channel_index: 0,
    };

    let results = run_multiwell_from_network(MultiWellFromNetworkInput { network, activation }, &micro)?;
    println!("Total outgoing rate: {}", results.total_outgoing_rate_constant);
    Ok(())

}
