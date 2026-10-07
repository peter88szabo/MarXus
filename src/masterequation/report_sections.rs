//! Sections of the human-readable master-equation report: the energetics of the deck (kcal/mol), and
//! the results of the steady-state and CSE solutions as named quantities at every (T, p), written in
//! the three views of `report_tables.rs` (by temperature, by pressure, temperature-pressure tables).
//!
//! Yields are given in percent, for every channel out of the network:
//! - steady state: the yields of the intermediate steady state (chemically activated, prompt), the
//!   thermal fate of the molecules thermalized in each well (final steady state with the Boltzmann
//!   distribution of that well as the source), the yields of the final steady state (both together), and
//!   their decomposition: prompt + stabilization x thermal fate, compared with the final steady state;
//! - CSE: the prompt branching of the reactant, the thermal fate of each well (absorbing chain of the well
//!   rate coefficients) and the long-time yields, direct + through the wells
//!   (`chemically_significant_eigenvalues::reactant_yields`).

use std::collections::HashMap;
use std::io::Write;

use crate::constants::CM1_TO_KCAL;

use super::chemical_activation_driver::{
    ConditionResult, PhenomenologicalConditionResult, ThermalConditionResult,
    WellFateConditionResult,
};
use super::chemical_activation_network::{
    ChannelDestination, ChemicalActivationNetwork, Conditions,
};
use super::chemical_activation_observables::ChemicalActivationResult;
use super::chemically_significant_eigenvalues::{reactant_yields, ReactantRates, ReactantYields};
use super::direct_time_integration::{TimeEvolution, TimePoint};
use super::mess_input::{MessBarrier, MessBarrierCore, MessDeck, TunnelingSpecification};
use super::report_tables::{
    sci, write_labelled_table, write_tables_by_pressure, write_tables_by_temperature,
    write_temperature_pressure_tables, Quantity, QuantityGroup,
};

// ------------------------------------------------------------------------------------------------------
// Energetics
// ------------------------------------------------------------------------------------------------------

/// Energy of a barrier in cm-1: the ZeroEnergy of its RRHO block in the deck (also for barrierless channels,
/// as the MESS output lists it).
fn barrier_energy_cm1(barrier: &MessBarrier) -> f64 {
    barrier.rrho.zero_energy_cm1
}

fn barrier_model(barrier: &MessBarrier) -> String {
    let mut model = if barrier.inverse_laplace_transform.is_some() {
        "barrierless, inverse Laplace transform".to_string()
    } else {
        match &barrier.core {
            MessBarrierCore::PhaseSpaceTheory { tst_level, .. } => {
                format!("barrierless, phase space theory ({tst_level:?})")
            }
            _ => "rigid transition state".to_string(),
        }
    };
    match &barrier.tunneling {
        Some(TunnelingSpecification::Eckart { .. }) => model.push_str(", Eckart tunneling"),
        Some(TunnelingSpecification::Unsupported { model: m }) => {
            model.push_str(&format!(", tunneling {m}"))
        }
        None => {}
    }
    model
}

/// Energetics of the deck in kcal/mol, relative to the Reactant if the deck names one: wells (ground
/// energy G, lowest barrier D out of the well), bimolecular species (ground energy) and barriers (ZeroEnergy
/// of the barrier in the deck).
pub fn write_energetics<W: Write>(out: &mut W, deck: &MessDeck) -> std::io::Result<()> {
    let reference = deck
        .global
        .reactant_name
        .as_ref()
        .and_then(|r| deck.bimolecular.get(r).map(|b| (r, b.ground_energy_cm1)));
    let zero = reference.map_or(0.0, |(_, e)| e);
    let kcal = |e_cm1: f64| format!("{:.2}", (e_cm1 - zero) * CM1_TO_KCAL);
    match reference {
        Some((name, _)) => writeln!(out, "Energetics (kcal/mol, relative to {name}):\n")?,
        None => writeln!(out, "Energetics (kcal/mol, energy zero of the deck):\n")?,
    }
    writeln!(
        out,
        "Wells (G - ground energy, D - lowest barrier out of the well):"
    )?;
    writeln!(out, "{:>10} {:>10} {:>10}", "Name", "G", "D")?;
    for name in &deck.well_order {
        let lowest = deck
            .barriers
            .iter()
            .filter(|b| b.left == *name || b.right == *name)
            .map(|b| barrier_energy_cm1(b))
            .fold(f64::INFINITY, f64::min);
        let d = if lowest.is_finite() {
            kcal(lowest)
        } else {
            "***".into()
        };
        writeln!(
            out,
            "{:>10} {:>10} {:>10}",
            name,
            kcal(deck.wells[name].zero_energy_cm1),
            d
        )?;
    }
    writeln!(out, "\nBimolecular (G - ground energy):")?;
    writeln!(out, "{:>10} {:>10}", "Name", "G")?;
    let mut bimolecular: Vec<_> = deck.bimolecular.values().collect();
    bimolecular.sort_by(|a, b| {
        b.ground_energy_cm1
            .total_cmp(&a.ground_energy_cm1)
            .then(a.name.cmp(&b.name))
    });
    for b in bimolecular {
        writeln!(out, "{:>10} {:>10}", b.name, kcal(b.ground_energy_cm1))?;
    }
    writeln!(
        out,
        "\nBarriers (H - ZeroEnergy of the barrier in the deck):"
    )?;
    writeln!(
        out,
        "{:>10} {:>10} {:>8} {:>8}   {}",
        "Name", "H", "From", "To", "Model"
    )?;
    for b in &deck.barriers {
        writeln!(
            out,
            "{:>10} {:>10} {:>8} {:>8}   {}",
            b.name,
            kcal(barrier_energy_cm1(b)),
            b.left,
            b.right,
            barrier_model(b)
        )?;
    }
    writeln!(out)
}

/// The chemical network of the run: every well (grains, collision parameters, sink) and every channel
/// (source and target, rate model from the deck: rigid transition state, phase space theory with its TST
/// level, inverse Laplace transform, tunneling), with the entrance channels marked.
pub fn write_network_summary<W: Write>(
    out: &mut W,
    deck: &MessDeck,
    network: &ChemicalActivationNetwork,
    entrance_channels: &[(usize, usize)],
) -> std::io::Result<()> {
    writeln!(out, "Wells (collisions: <dE_down>(T) = <dE_down>(T_ref) (T/T_ref)^n; Lennard-Jones collision frequency):")?;
    writeln!(
        out,
        "{:>10} {:>8} {:>18}   {}",
        "Name", "grains", "absolute grains", "collisions; sink"
    )?;
    for well in &network.wells {
        let et = &well.energy_transfer;
        let lj = &well.lennard_jones;
        let range = format!(
            "{} .. {}",
            well.bottom_offset_grains,
            well.bottom_offset_grains + well.grain_count() as isize - 1
        );
        let sink = if well.bimolecular_sink_s_inv > 0.0 {
            format!(
                "; escape (bimolecular sink) {} 1/s",
                sci(well.bimolecular_sink_s_inv)
            )
        } else {
            String::new()
        };
        writeln!(
            out,
            "{:>10} {:>8} {:>18}   <dE_down> = {} cm-1 (T/{} K)^{}; sigma {:.3} A, epsilon {:.2} K, reduced mass {:.3} amu{sink}",
            well.name,
            well.grain_count(),
            range,
            et.mean_down_at_reference_cm1,
            et.reference_temperature_kelvin,
            et.temperature_exponent,
            lj.sigma_angstrom,
            lj.epsilon_kelvin,
            lj.reduced_mass_amu,
        )?;
    }
    writeln!(
        out,
        "\nChannels (k(E) of every channel; rate model from the deck):"
    )?;
    writeln!(
        out,
        "{:>10} {:>8} {:>8}   {}",
        "Name", "From", "To", "Model"
    )?;
    for (w, well) in network.wells.iter().enumerate() {
        for (c, channel) in well.channels.iter().enumerate() {
            let to = match &channel.destination {
                ChannelDestination::Products { name } => name.clone(),
                ChannelDestination::Well { index } => network.wells[*index].name.clone(),
            };
            let model = deck
                .barriers
                .iter()
                .find(|b| b.name == channel.name)
                .map_or_else(|| "(not a deck barrier)".to_string(), barrier_model);
            let entrance = if entrance_channels.contains(&(w, c)) {
                "; entrance of the reactant (source)"
            } else {
                ""
            };
            writeln!(
                out,
                "{:>10} {:>8} {:>8}   {model}{entrance}",
                channel.name, well.name, to
            )?;
        }
    }
    writeln!(out)
}

// ------------------------------------------------------------------------------------------------------
// Quantities at every condition
// ------------------------------------------------------------------------------------------------------

fn key(t: f64, p: f64) -> (u64, u64) {
    (t.to_bits(), p.to_bits())
}

fn condition_key(c: &Conditions) -> (u64, u64) {
    key(c.temperature_kelvin, c.pressure_torr)
}

/// A quantity from the results indexed by condition; None where a condition has no result.
fn quantity<R>(
    name: impl Into<String>,
    temperatures: &[f64],
    pressures: &[f64],
    index: &HashMap<(u64, u64), &R>,
    value: impl Fn(&R, f64) -> Option<f64>,
) -> Quantity {
    Quantity {
        name: name.into(),
        values: temperatures
            .iter()
            .map(|&t| {
                pressures
                    .iter()
                    .map(|&p| index.get(&key(t, p)).and_then(|r| value(r, t)))
                    .collect()
            })
            .collect(),
    }
}

/// Exits of the network: a product channel of a well, the stabilization of a well (intermediate steady
/// state), or the bimolecular sink of a well.
#[derive(Debug, Clone, Copy)]
enum Exit {
    Product { well: usize, channel: usize },
    Stabilization(usize),
    Sink(usize),
}

/// Name of a channel: `W->X`, X the product or the target well.
fn channel_name(network: &ChemicalActivationNetwork, well: usize, channel: usize) -> String {
    let w = &network.wells[well];
    let target = match &w.channels[channel].destination {
        ChannelDestination::Products { name } => name.clone(),
        ChannelDestination::Well { index } => network.wells[*index].name.clone(),
    };
    format!("{}->{target}", w.name)
}

/// Unique names: a name that occurs more than once gets the channel name appended.
fn unique_channel_names(
    network: &ChemicalActivationNetwork,
    channels: &[(usize, usize)],
) -> Vec<String> {
    let names: Vec<String> = channels
        .iter()
        .map(|&(w, c)| channel_name(network, w, c))
        .collect();
    names
        .iter()
        .zip(channels)
        .map(|(n, &(w, c))| {
            if names.iter().filter(|m| *m == n).count() > 1 {
                format!("{n}[{}]", network.wells[w].channels[c].name)
            } else {
                n.clone()
            }
        })
        .collect()
}

fn all_channels(network: &ChemicalActivationNetwork) -> Vec<(usize, usize)> {
    network
        .wells
        .iter()
        .enumerate()
        .flat_map(|(w, well)| (0..well.channels.len()).map(move |c| (w, c)))
        .collect()
}

fn product_channels(network: &ChemicalActivationNetwork) -> Vec<(usize, usize)> {
    all_channels(network)
        .into_iter()
        .filter(|&(w, c)| {
            matches!(
                network.wells[w].channels[c].destination,
                ChannelDestination::Products { .. }
            )
        })
        .collect()
}

fn product_name(network: &ChemicalActivationNetwork, well: usize, channel: usize) -> Option<&str> {
    match &network.wells[well].channels[channel].destination {
        ChannelDestination::Products { name } => Some(name),
        ChannelDestination::Well { .. } => None,
    }
}

/// Product channels and sinks, with their names; with `stabilization` also the stabilization of every well.
fn exits(network: &ChemicalActivationNetwork, stabilization: bool) -> Vec<(String, Exit)> {
    let products = product_channels(network);
    let mut list: Vec<(String, Exit)> = unique_channel_names(network, &products)
        .into_iter()
        .zip(&products)
        .map(|(n, &(well, channel))| (n, Exit::Product { well, channel }))
        .collect();
    if stabilization {
        list.extend(
            network
                .wells
                .iter()
                .enumerate()
                .map(|(w, well)| (format!("stab({})", well.name), Exit::Stabilization(w))),
        );
    }
    list.extend(
        network
            .wells
            .iter()
            .enumerate()
            .filter(|(_, w)| w.bimolecular_sink_s_inv > 0.0)
            .map(|(w, well)| (format!("escape({})", well.name), Exit::Sink(w))),
    );
    list
}

/// Targets of the reactant `r`: its bimolecular products `r->P` (all product channels to P summed; the
/// return to `r` excluded), the sinks `r->escape(W)`, and the wells `r->W` (stabilization), in this order.
fn reactant_targets(
    network: &ChemicalActivationNetwork,
    exits: &[(String, Exit)],
    r: &str,
) -> Vec<(String, Vec<Exit>)> {
    let mut targets: Vec<(String, Vec<Exit>)> = Vec::new();
    for (_, e) in exits {
        let name = match e {
            Exit::Product { well, channel } => match product_name(network, *well, *channel) {
                Some(p) if p == r => continue,
                Some(p) => format!("{r}->{p}"),
                None => continue,
            },
            Exit::Sink(w) => format!("{r}->escape({})", network.wells[*w].name),
            Exit::Stabilization(_) => continue,
        };
        match targets.iter_mut().find(|(n, _)| *n == name) {
            Some((_, members)) => members.push(*e),
            None => targets.push((name, vec![*e])),
        }
    }
    for (_, e) in exits {
        if let Exit::Stabilization(w) = e {
            targets.push((format!("{r}->{}", network.wells[*w].name), vec![*e]));
        }
    }
    targets
}

fn exit_value(result: &ChemicalActivationResult, exit: Exit) -> f64 {
    match exit {
        Exit::Product { well, channel } => result
            .channels
            .iter()
            .find(|c| c.well == well && c.channel == channel)
            .map_or(f64::NAN, |c| c.flux),
        Exit::Stabilization(w) => result.wells[w].stabilization_yield,
        Exit::Sink(w) => result.wells[w].bimolecular_sink_yield,
    }
}

fn ca_rate(result: &ChemicalActivationResult, well: usize, channel: usize) -> f64 {
    result
        .channels
        .iter()
        .find(|c| c.well == well && c.channel == channel)
        .map_or(f64::NAN, |c| c.ca_rate_constant_s_inv)
}

// ------------------------------------------------------------------------------------------------------
// Steady state
// ------------------------------------------------------------------------------------------------------

/// Groups of one steady-state solution: yields (% of the formed adducts), yields without the return to
/// the reactant (% of the net reaction), chemical-activation rate coefficients k_ca (1/s) and, with
/// `k_inf`, the bimolecular rate coefficients k(R -> X) = k_inf Phi_X (cm3/s).
pub fn steady_state_groups(
    network: &ChemicalActivationNetwork,
    temperatures: &[f64],
    pressures: &[f64],
    results: &[ConditionResult],
    reactant: Option<&str>,
    k_inf: Option<&dyn Fn(f64) -> f64>,
) -> Vec<QuantityGroup> {
    let index: HashMap<_, _> = results
        .iter()
        .map(|r| (condition_key(&r.conditions), r))
        .collect();
    let stabilization = results
        .iter()
        .any(|r| r.result.wells.iter().any(|w| w.stabilization_yield != 0.0));
    let exits = exits(network, stabilization);
    let mut groups = Vec::new();

    groups.push(QuantityGroup {
        title: "Yields (% of the formed adducts)".into(),
        quantities: exits
            .iter()
            .map(|(n, e)| {
                quantity(
                    n.clone(),
                    temperatures,
                    pressures,
                    &index,
                    |r: &ConditionResult, _| Some(100.0 * exit_value(&r.result, *e)),
                )
            })
            .collect(),
    });

    if let Some(r_name) = reactant {
        let is_back = |e: &Exit| matches!(e, Exit::Product { well, channel } if product_name(network, *well, *channel) == Some(r_name));
        let back: Vec<Exit> = exits
            .iter()
            .map(|(_, e)| *e)
            .filter(|e| is_back(e))
            .collect();
        if !back.is_empty() {
            let back_yield =
                |r: &ConditionResult| back.iter().map(|e| exit_value(&r.result, *e)).sum::<f64>();
            groups.push(QuantityGroup {
                title: format!("Yields without the return to {r_name} (% of the net reaction)"),
                quantities: exits
                    .iter()
                    .filter(|(_, e)| !is_back(e))
                    .map(|(n, e)| {
                        quantity(
                            n.clone(),
                            temperatures,
                            pressures,
                            &index,
                            |r: &ConditionResult, _| {
                                Some(100.0 * exit_value(&r.result, *e) / (1.0 - back_yield(r)))
                            },
                        )
                    })
                    .collect(),
            });
        }
    }

    let channels = all_channels(network);
    groups.push(QuantityGroup {
        title: "Chemical-activation rate coefficients k_ca (1/s)".into(),
        quantities: unique_channel_names(network, &channels)
            .into_iter()
            .zip(&channels)
            .map(|(n, &(w, c))| {
                quantity(
                    n,
                    temperatures,
                    pressures,
                    &index,
                    |r: &ConditionResult, _| Some(ca_rate(&r.result, w, c)),
                )
            })
            .collect(),
    });

    if let (Some(r_name), Some(k_inf)) = (reactant, k_inf) {
        // Bimolecular-to-bimolecular (chemically activated products, summed over the channels that form them,
        // and the sinks) and bimolecular-to-well (stabilization) rate coefficients k(R -> X) = k_inf Phi_X
        // (PR03 eq. 44) and yields. In the final steady state there is no stabilization, and the products
        // include the thermal reaction of the stabilized adducts: "overall".
        let targets = reactant_targets(network, &exits, r_name);
        let overall = !stabilization;
        let back = |r: &ConditionResult| {
            exits
                .iter()
                .filter(|(_, e)| matches!(e, Exit::Product { well, channel } if product_name(network, *well, *channel) == Some(r_name)))
                .map(|(_, e)| exit_value(&r.result, *e))
                .sum::<f64>()
        };
        let target_value = |r: &ConditionResult, members: &[Exit]| {
            members
                .iter()
                .map(|e| exit_value(&r.result, *e))
                .sum::<f64>()
        };
        let (bb_title, bb_yield_title) = if overall {
            (
                format!("Bimolecular-to-bimolecular rate coefficients, overall (chemical activation + thermal reaction of the stabilized adducts): k({r_name} -> P) = k_inf Phi_P (cm^3/s)"),
                format!("Bimolecular-to-bimolecular yields, overall (chemical activation + thermal) (% of the net reaction of {r_name})"),
            )
        } else {
            (
                format!("Bimolecular-to-bimolecular rate coefficients (chemical activation): k({r_name} -> P) = k_inf Phi_P (cm^3/s)"),
                format!("Bimolecular-to-bimolecular yields (chemical activation) (% of the net reaction of {r_name})"),
            )
        };
        let bimolecular: Vec<&(String, Vec<Exit>)> = targets
            .iter()
            .filter(|(_, m)| !matches!(m[0], Exit::Stabilization(_)))
            .collect();
        let wells: Vec<&(String, Vec<Exit>)> = targets
            .iter()
            .filter(|(_, m)| matches!(m[0], Exit::Stabilization(_)))
            .collect();
        let rate_group = |title: String, list: &[&(String, Vec<Exit>)]| QuantityGroup {
            title,
            quantities: list
                .iter()
                .map(|(n, m)| {
                    quantity(
                        n.clone(),
                        temperatures,
                        pressures,
                        &index,
                        |r: &ConditionResult, t| Some(k_inf(t) * target_value(r, m)),
                    )
                })
                .collect(),
        };
        let yield_group = |title: String, list: &[&(String, Vec<Exit>)]| QuantityGroup {
            title,
            quantities: list
                .iter()
                .map(|(n, m)| {
                    quantity(
                        n.clone(),
                        temperatures,
                        pressures,
                        &index,
                        |r: &ConditionResult, _| Some(100.0 * target_value(r, m) / (1.0 - back(r))),
                    )
                })
                .collect(),
        };
        if !bimolecular.is_empty() {
            groups.push(rate_group(bb_title, &bimolecular));
            groups.push(yield_group(bb_yield_title, &bimolecular));
        }
        if !wells.is_empty() {
            groups.push(rate_group(
                format!("Bimolecular-to-well rate coefficients (stabilization): k({r_name} -> W) = k_inf Phi_stab,W (cm^3/s)"),
                &wells,
            ));
            groups.push(yield_group(format!("Bimolecular-to-well yields (stabilization) (% of the net reaction of {r_name})"), &wells));
        }
        groups.push(QuantityGroup {
            title: format!("Capture, return and net reaction of {r_name}: k_inf, k_inf Phi_return, k_inf (1 - Phi_return) (cm^3/s)"),
            quantities: vec![
                quantity("capture", temperatures, pressures, &index, |_: &ConditionResult, t| Some(k_inf(t))),
                quantity("return", temperatures, pressures, &index, |r: &ConditionResult, t| Some(k_inf(t) * back(r))),
                quantity("net", temperatures, pressures, &index, |r: &ConditionResult, t| Some(k_inf(t) * (1.0 - back(r)))),
            ],
        });
    }
    groups
}

/// One group per well: the thermal fate (%) of the molecules thermalized in that well.
pub fn well_fate_groups(
    network: &ChemicalActivationNetwork,
    temperatures: &[f64],
    pressures: &[f64],
    fates: &[WellFateConditionResult],
) -> Vec<QuantityGroup> {
    let exits = exits(network, false);
    network
        .wells
        .iter()
        .enumerate()
        .map(|(w, well)| {
            let index: HashMap<_, _> = fates
                .iter()
                .filter(|f| f.well == w)
                .map(|f| (condition_key(&f.conditions), f))
                .collect();
            QuantityGroup {
                title: format!(
                    "Thermal fate of the molecules thermalized in {} (%)",
                    well.name
                ),
                quantities: exits
                    .iter()
                    .map(|(n, e)| {
                        quantity(
                            n.clone(),
                            temperatures,
                            pressures,
                            &index,
                            |f: &WellFateConditionResult, _| {
                                Some(100.0 * exit_value(&f.result, *e))
                            },
                        )
                    })
                    .collect(),
            }
        })
        .collect()
}

/// Chemical activation and thermal reaction separately and together (% of the formed adducts): the prompt
/// yields of the intermediate steady state, the yields through the stabilized wells (stabilization yield
/// x thermal fate of the well), their sum, the final steady state, and the difference of the two.
pub fn chemical_activation_and_thermal_groups(
    network: &ChemicalActivationNetwork,
    temperatures: &[f64],
    pressures: &[f64],
    intermediate: &[ConditionResult],
    fates: &[WellFateConditionResult],
    final_steady_state: &[ConditionResult],
) -> Vec<QuantityGroup> {
    let exits = exits(network, false);
    let prompt: HashMap<_, _> = intermediate
        .iter()
        .map(|r| (condition_key(&r.conditions), r))
        .collect();
    let total: HashMap<_, _> = final_steady_state
        .iter()
        .map(|r| (condition_key(&r.conditions), r))
        .collect();
    let fate: HashMap<_, _> = fates
        .iter()
        .map(|f| ((condition_key(&f.conditions), f.well), f))
        .collect();
    // Values of one exit at one condition: (prompt, through the wells, final), None if any part is missing.
    let parts = |t: f64, p: f64, e: Exit| -> Option<(f64, f64, f64)> {
        let k = key(t, p);
        let r = prompt.get(&k)?;
        let mut via = 0.0;
        for (w, well) in r.result.wells.iter().enumerate() {
            via += well.stabilization_yield * exit_value(&fate.get(&(k, w))?.result, e);
        }
        Some((
            exit_value(&r.result, e),
            via,
            exit_value(&total.get(&k)?.result, e),
        ))
    };
    let make = |title: &str, f: &dyn Fn((f64, f64, f64)) -> f64| QuantityGroup {
        title: title.into(),
        quantities: exits
            .iter()
            .map(|(n, e)| Quantity {
                name: n.clone(),
                values: temperatures
                    .iter()
                    .map(|&t| pressures.iter().map(|&p| parts(t, p, *e).map(f)).collect())
                    .collect(),
            })
            .collect(),
    };
    vec![
        make("Chemically activated (prompt) yields: intermediate steady state (% of the formed adducts)", &|(a, _, _)| 100.0 * a),
        make("Through the stabilized wells: stabilization yield x thermal fate of the well (% of the formed adducts)", &|(_, b, _)| 100.0 * b),
        make("Prompt + through the stabilized wells (% of the formed adducts)", &|(a, b, _)| 100.0 * (a + b)),
        make("All together: final steady state (% of the formed adducts)", &|(_, _, c)| 100.0 * c),
        make("Difference: prompt + through the stabilized wells - final steady state (percentage points)", &|(a, b, c)| 100.0 * (a + b - c)),
    ]
}

// ------------------------------------------------------------------------------------------------------
// Thermal rate coefficients of the final steady state
// ------------------------------------------------------------------------------------------------------

/// Groups of the thermal rate coefficients of the final steady state: k_uni, lambda_1 and the channel
/// rate coefficients (1/s), the thermal branching (%), the high-pressure rate coefficients (1/s) and the
/// diagnostics of the eigenpair.
pub fn thermal_groups(
    network: &ChemicalActivationNetwork,
    temperatures: &[f64],
    pressures: &[f64],
    results: &[ThermalConditionResult],
) -> Vec<QuantityGroup> {
    let index: HashMap<_, _> = results
        .iter()
        .map(|r| (condition_key(&r.conditions), r))
        .collect();
    let channels = all_channels(network);
    let names = unique_channel_names(network, &channels);
    let channel_rate = |r: &ThermalConditionResult, w: usize, c: usize, high_pressure: bool| {
        r.thermal
            .channels
            .iter()
            .find(|x| x.well == w && x.channel == c)
            .map(|x| {
                if high_pressure {
                    x.high_pressure_rate_s_inv
                } else {
                    x.thermal_rate_s_inv
                }
            })
    };
    let sinks: Vec<usize> = (0..network.wells.len())
        .filter(|&w| network.wells[w].bimolecular_sink_s_inv > 0.0)
        .collect();
    let sink_name = |w: usize| format!("escape({})", network.wells[w].name);

    let mut rates = vec![
        quantity(
            "k_uni",
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| Some(r.thermal.k_uni_s_inv),
        ),
        quantity(
            "lambda_1",
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| Some(r.thermal.lambda_1_s_inv),
        ),
    ];
    for (n, &(w, c)) in names.iter().zip(&channels) {
        rates.push(quantity(
            n.clone(),
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| channel_rate(r, w, c, false),
        ));
    }
    for &w in &sinks {
        rates.push(quantity(
            sink_name(w),
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| Some(r.thermal.sink_rates_s_inv[w]),
        ));
    }

    let mut branching = Vec::new();
    for (n, &(w, c)) in names.iter().zip(&channels) {
        if product_name(network, w, c).is_some() {
            branching.push(quantity(
                n.clone(),
                temperatures,
                pressures,
                &index,
                |r: &ThermalConditionResult, _| {
                    channel_rate(r, w, c, false).map(|k| 100.0 * k / r.thermal.k_uni_s_inv)
                },
            ));
        }
    }
    for &w in &sinks {
        branching.push(quantity(
            sink_name(w),
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| {
                Some(100.0 * r.thermal.sink_rates_s_inv[w] / r.thermal.k_uni_s_inv)
            },
        ));
    }

    let high_pressure = names
        .iter()
        .zip(&channels)
        .map(|(n, &(w, c))| {
            quantity(
                n.clone(),
                temperatures,
                pressures,
                &index,
                |r: &ThermalConditionResult, _| channel_rate(r, w, c, true),
            )
        })
        .collect();

    let diagnostics = vec![
        quantity(
            "sum rule",
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| Some(r.thermal.sum_rule_relative_deviation),
        ),
        quantity(
            "lambda_2/k_uni",
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| Some(r.thermal.lambda_2_s_inv / r.thermal.k_uni_s_inv),
        ),
        quantity(
            "floor (1/s)",
            temperatures,
            pressures,
            &index,
            |r: &ThermalConditionResult, _| Some(r.thermal.precision_floor_s_inv),
        ),
    ];

    vec![
        QuantityGroup { title: "Thermal rate coefficients (1/s)".into(), quantities: rates },
        QuantityGroup { title: "Thermal branching (% of the thermal decay k_uni)".into(), quantities: branching },
        QuantityGroup { title: "High-pressure rate coefficients k_inf (1/s)".into(), quantities: high_pressure },
        QuantityGroup {
            title: "Diagnostics of the thermal eigenpair (sum rule |lambda_1 - k_uni|/k_uni; lambda_2/k_uni; precision floor eps max S_ii)".into(),
            quantities: diagnostics,
        },
    ]
}

// ------------------------------------------------------------------------------------------------------
// CSE
// ------------------------------------------------------------------------------------------------------

/// Groups of the CSE solution: the rate coefficients from every species (1/s) and from the reactant (cm3/s),
/// the yields of the reactant and the thermal fates of the species (%), and the diagnostics. Where wells are
/// merged (Georgievskii et al. 2013, Sec. IV) the species differ between conditions: the rows are all species
/// of all conditions, by name, and empty where a species does not exist.
pub fn cse_groups(
    temperatures: &[f64],
    pressures: &[f64],
    results: &[PhenomenologicalConditionResult],
) -> Vec<QuantityGroup> {
    let Some(first) = results.first() else {
        return Vec::new();
    };
    let mut wells: Vec<String> = Vec::new();
    for r in results {
        for w in &r.rates.wells {
            if !wells.contains(w) {
                wells.push(w.clone());
            }
        }
    }
    // Index of a species at one condition.
    let at = |r: &PhenomenologicalConditionResult, name: &str| r.rates.wells.iter().position(|w| w == name);
    let bimolecular = first.rates.bimolecular.clone();
    let index: HashMap<_, _> = results
        .iter()
        .map(|r| (condition_key(&r.conditions), r))
        .collect();
    let mut groups = Vec::new();

    for well in &wells {
        let mut quantities = Vec::new();
        for target in wells.iter().filter(|&t| t != well) {
            // Only species that coexist at some condition.
            if !results.iter().any(|r| at(r, well).is_some() && at(r, target).is_some()) {
                continue;
            }
            quantities.push(quantity(
                format!("{well}->{target}"),
                temperatures,
                pressures,
                &index,
                |r: &PhenomenologicalConditionResult, _| {
                    Some(r.rates.well_to_well_s_inv[at(r, well)?][at(r, target)?])
                },
            ));
        }
        for (nu, target) in bimolecular.iter().enumerate() {
            quantities.push(quantity(
                format!("{well}->{target}"),
                temperatures,
                pressures,
                &index,
                |r: &PhenomenologicalConditionResult, _| {
                    Some(r.rates.well_to_bimolecular_s_inv[at(r, well)?][nu])
                },
            ));
        }
        quantities.push(quantity(
            format!("{well} loss"),
            temperatures,
            pressures,
            &index,
            |r: &PhenomenologicalConditionResult, _| {
                let i = at(r, well)?;
                Some(r.rates.well_to_well_s_inv[i][i])
            },
        ));
        groups.push(QuantityGroup {
            title: format!("Rate coefficients from {well} (1/s)"),
            quantities,
        });
    }

    if let Some(reactant) = &first.rates.reactant {
        let r_name = reactant.name.clone();
        let r_index = bimolecular.iter().position(|b| *b == r_name);
        let rate = |name: String, pick: Box<dyn Fn(&PhenomenologicalConditionResult, &ReactantRates) -> Option<f64>>| {
            quantity(
                name,
                temperatures,
                pressures,
                &index,
                move |r: &PhenomenologicalConditionResult, _| {
                    r.rates.reactant.as_ref().and_then(|x| pick(r, x))
                },
            )
        };
        // Bimolecular-to-bimolecular (chemical activation, G13 eq. 21) and bimolecular-to-well (stabilization,
        // G13 eq. 28) rate coefficients of the reactant.
        let to_bimolecular: Vec<Quantity> = bimolecular
            .iter()
            .enumerate()
            .filter(|&(nu, _)| Some(nu) != r_index)
            .map(|(nu, target)| {
                rate(
                    format!("{r_name}->{target}"),
                    Box::new(move |_, x: &ReactantRates| Some(x.to_bimolecular_cm3_s[nu])),
                )
            })
            .collect();
        if !to_bimolecular.is_empty() {
            groups.push(QuantityGroup {
                title: format!("Bimolecular-to-bimolecular rate coefficients (chemical activation; G13 eq. 21): k({r_name} -> P) (cm^3/s)"),
                quantities: to_bimolecular,
            });
        }
        groups.push(QuantityGroup {
            title: format!("Bimolecular-to-well rate coefficients (stabilization; G13 eq. 28): k({r_name} -> W) (cm^3/s)"),
            quantities: wells
                .iter()
                .map(|well| {
                    let name = well.clone();
                    rate(
                        format!("{r_name}->{well}"),
                        Box::new(move |r, x: &ReactantRates| Some(x.to_well_cm3_s[at(r, &name)?])),
                    )
                })
                .collect(),
        });
        if let Some(ri) = r_index {
            groups.push(QuantityGroup {
                title: format!("Capture, return and net reaction of {r_name} (cm^3/s)"),
                quantities: vec![
                    rate(
                        "capture".into(),
                        Box::new(|_, x: &ReactantRates| Some(x.capture_cm3_s)),
                    ),
                    rate(
                        "return".into(),
                        Box::new(move |_, x: &ReactantRates| Some(x.to_bimolecular_cm3_s[ri])),
                    ),
                    rate(
                        "net".into(),
                        Box::new(move |_, x: &ReactantRates| {
                            Some(x.capture_cm3_s - x.to_bimolecular_cm3_s[ri])
                        }),
                    ),
                ],
            });
        }

        // Yields of the reactant, computed once per condition.
        let yields: HashMap<(u64, u64), (ReactantYields, &PhenomenologicalConditionResult)> = results
            .iter()
            .filter_map(|r| {
                reactant_yields(&r.rates)
                    .and_then(|y| y.ok())
                    .map(|y| (condition_key(&r.conditions), (y, r)))
            })
            .collect();
        let yield_group = |title: String,
                           names: Vec<String>,
                           pick: &dyn Fn(&ReactantYields, &PhenomenologicalConditionResult, usize) -> Option<f64>| {
            QuantityGroup {
                title,
                quantities: names
                    .into_iter()
                    .enumerate()
                    .map(|(q, name)| Quantity {
                        name,
                        values: temperatures
                            .iter()
                            .map(|&t| {
                                pressures
                                    .iter()
                                    .map(|&p| yields.get(&key(t, p)).and_then(|(y, r)| pick(y, r, q)).map(|x| 100.0 * x))
                                    .collect()
                            })
                            .collect(),
                    })
                    .collect(),
            }
        };
        let others: Vec<String> = bimolecular
            .iter()
            .filter(|b| **b != r_name)
            .cloned()
            .collect();
        // The prompt branching of the net reaction, split into its bimolecular and its well part (the species of
        // the condition first).
        if !others.is_empty() {
            groups.push(yield_group(
                format!("Bimolecular-to-bimolecular yields (chemical activation) (% of the net reaction of {r_name})"),
                others.iter().map(|x| format!("{r_name}->{x}")).collect(),
                &|y, r, q| Some(y.prompt_branching[r.rates.wells.len() + q]),
            ));
        }
        groups.push(yield_group(
            format!(
                "Bimolecular-to-well yields (stabilization) (% of the net reaction of {r_name})"
            ),
            wells.iter().map(|x| format!("{r_name}->{x}")).collect(),
            &|y, r, q| Some(y.prompt_branching[at(r, &wells[q])?]),
        ));
        for well in &wells {
            let names = bimolecular.iter().map(|x| format!("{well}->{x}")).collect();
            groups.push(yield_group(
                format!("Thermal fate of {well} (%)"),
                names,
                &|y, r, q| Some(y.well_fates[at(r, well)?][q]),
            ));
        }
        let channels: Vec<String> = others.iter().map(|x| format!("{r_name}->{x}")).collect();
        groups.push(yield_group(
            format!("Long-time yields, direct (chemically activated) (% of the eventual net reaction of {r_name})"),
            channels.clone(),
            &|y, _, q| Some(y.direct[q]),
        ));
        groups.push(yield_group(
            format!("Long-time yields, through the wells (thermal) (% of the eventual net reaction of {r_name})"),
            channels.clone(),
            &|y, _, q| Some(y.via_wells[q]),
        ));
        groups.push(yield_group(
            format!("Long-time yields, total (% of the eventual net reaction of {r_name})"),
            channels,
            &|y, _, q| Some(y.total[q]),
        ));
    }

    groups.push(QuantityGroup {
        title: "Diagnostics of the CSE solution (species: number of kinetic species, fewer than the wells where wells \
                are merged; separation Lambda_N/Lambda_N+1; loss balance; detailed balance)"
            .into(),
        quantities: vec![
            quantity("species", temperatures, pressures, &index, |r: &PhenomenologicalConditionResult, _| {
                Some(r.rates.wells.len() as f64)
            }),
            quantity("separation", temperatures, pressures, &index, |r: &PhenomenologicalConditionResult, _| {
                r.rates.chemical_eigenvalues_s_inv.last().map(|l| l / r.rates.relaxation_eigenvalue_s_inv)
            }),
            quantity("loss balance", temperatures, pressures, &index, |r: &PhenomenologicalConditionResult, _| {
                Some(r.rates.loss_balance_max_deviation)
            }),
            quantity("detailed bal.", temperatures, pressures, &index, |r: &PhenomenologicalConditionResult, _| {
                Some(r.rates.detailed_balance_max_deviation)
            }),
        ],
    });
    groups
}

/// Species-to-species tables of the CSE solution, one per condition, with the eigenvalues, the capture
/// and the warnings above each table.
pub fn write_cse_species_tables<W: Write>(
    out: &mut W,
    results: &[PhenomenologicalConditionResult],
) -> std::io::Result<()> {
    for r in results {
        let x = &r.rates;
        writeln!(
            out,
            "Temperature = {} K    Pressure = {} torr\n",
            r.conditions.temperature_kelvin, r.conditions.pressure_torr
        )?;
        let eigenvalues: Vec<String> = x
            .chemical_eigenvalues_s_inv
            .iter()
            .map(|&l| sci(l))
            .collect();
        writeln!(
            out,
            "  chemical eigenvalues (1/s):          {}",
            eigenvalues.join("  ")
        )?;
        writeln!(
            out,
            "  lowest relaxation eigenvalue (1/s):  {}   separation Lambda_N/Lambda_N+1: {}",
            sci(x.relaxation_eigenvalue_s_inv),
            sci(x
                .chemical_eigenvalues_s_inv
                .last()
                .copied()
                .unwrap_or(f64::NAN)
                / x.relaxation_eigenvalue_s_inv)
        )?;
        let mut rows: Vec<(String, Vec<Option<f64>>)> = x
            .wells
            .iter()
            .enumerate()
            .map(|(i, w)| {
                (
                    w.clone(),
                    x.well_to_well_s_inv[i]
                        .iter()
                        .chain(&x.well_to_bimolecular_s_inv[i])
                        .map(|&k| Some(k))
                        .collect(),
                )
            })
            .collect();
        if let Some(reactant) = &x.reactant {
            let r_index = x.bimolecular.iter().position(|b| *b == reactant.name);
            if let Some(ri) = r_index {
                writeln!(
                    out,
                    "  {} (cm^3/s): capture {}   return {}   net reaction {}",
                    reactant.name,
                    sci(reactant.capture_cm3_s),
                    sci(reactant.to_bimolecular_cm3_s[ri]),
                    sci(reactant.capture_cm3_s - reactant.to_bimolecular_cm3_s[ri])
                )?;
            }
            let mut row: Vec<Option<f64>> =
                reactant.to_well_cm3_s.iter().map(|&k| Some(k)).collect();
            row.extend(
                reactant
                    .to_bimolecular_cm3_s
                    .iter()
                    .enumerate()
                    .map(|(nu, &k)| {
                        Some(if Some(nu) == r_index {
                            reactant.capture_cm3_s - k
                        } else {
                            k
                        })
                    }),
            );
            rows.push((reactant.name.clone(), row));
        }
        for warning in &x.warnings {
            writeln!(out, "  warning: {warning}")?;
        }
        writeln!(out)?;
        let columns: Vec<String> = x.wells.iter().chain(&x.bimolecular).cloned().collect();
        write_labelled_table(out, "From\\To", &columns, &rows)?;
    }
    Ok(())
}

// ------------------------------------------------------------------------------------------------------
// Direct time integration
// ------------------------------------------------------------------------------------------------------

/// The time evolution of every condition: a table with the output times as rows and the population of
/// every well and the yield of every exit (both in % of the formed adducts) as columns, with the work of
/// the integrator above it.
pub fn write_time_evolution_tables<W: Write>(
    out: &mut W,
    network: &ChemicalActivationNetwork,
    evolutions: &[TimeEvolution],
) -> std::io::Result<()> {
    for e in evolutions {
        writeln!(
            out,
            "Temperature = {} K    Pressure = {} torr\n",
            e.conditions.temperature_kelvin, e.conditions.pressure_torr
        )?;
        let st = &e.statistics;
        writeln!(
            out,
            "  integrator: {} steps ({} accepted, {} rejected), {} function evaluations, {} factorizations computed\n",
            st.steps, st.accepted_steps, st.rejected_steps, st.function_evaluations, e.factorizations_computed
        )?;
        let mut columns: Vec<String> = network
            .wells
            .iter()
            .map(|w| format!("N({})", w.name))
            .collect();
        columns.extend(e.exits.iter().cloned());
        let rows: Vec<(String, Vec<Option<f64>>)> = e
            .points
            .iter()
            .map(|p| {
                let values = p
                    .well_populations
                    .iter()
                    .chain(&p.exit_yields)
                    .map(|&v| Some(100.0 * v))
                    .collect();
                (sci(p.time_s), values)
            })
            .collect();
        write_labelled_table(out, "t(s)", &columns, &rows)?;
    }
    Ok(())
}

/// Groups of the time integration at its last output time: the exit yields (% of the formed adducts), the
/// populations left in the wells (%), with `reactant` the yields without the return to it (% of the net
/// reaction; equal to the final steady state for a completed pulse), and the conserved total (%).
pub fn time_integration_groups(
    network: &ChemicalActivationNetwork,
    temperatures: &[f64],
    pressures: &[f64],
    evolutions: &[TimeEvolution],
    reactant: Option<&str>,
    k_inf: Option<&dyn Fn(f64) -> f64>,
) -> Vec<QuantityGroup> {
    let Some(first) = evolutions.first() else {
        return Vec::new();
    };
    let index: HashMap<_, _> = evolutions
        .iter()
        .map(|e| (condition_key(&e.conditions), e))
        .collect();
    let last = |e: &TimeEvolution| e.points.last().cloned();
    let mut groups = vec![
        QuantityGroup {
            title: "Yields at the last output time (% of the formed adducts)".into(),
            quantities: first
                .exits
                .iter()
                .enumerate()
                .map(|(x, name)| {
                    quantity(
                        name.clone(),
                        temperatures,
                        pressures,
                        &index,
                        |e: &TimeEvolution, _| last(e).map(|p| 100.0 * p.exit_yields[x]),
                    )
                })
                .collect(),
        },
        QuantityGroup {
            title:
                "Populations left in the wells at the last output time (% of the formed adducts)"
                    .into(),
            quantities: network
                .wells
                .iter()
                .enumerate()
                .map(|(w, well)| {
                    quantity(
                        format!("N({})", well.name),
                        temperatures,
                        pressures,
                        &index,
                        |e: &TimeEvolution, _| last(e).map(|p| 100.0 * p.well_populations[w]),
                    )
                })
                .collect(),
        },
    ];
    if let Some(r) = reactant {
        let back: Vec<usize> = first
            .exits
            .iter()
            .enumerate()
            .filter(|(_, n)| n.ends_with(&format!("->{r}")))
            .map(|(x, _)| x)
            .collect();
        if !back.is_empty() {
            let net = |p: &TimePoint| {
                p.exit_yields
                    .iter()
                    .enumerate()
                    .filter(|(x, _)| !back.contains(x))
                    .map(|(_, y)| y)
                    .sum::<f64>()
            };
            groups.push(QuantityGroup {
                title: format!("Yields without the return to {r} at the last output time (% of the net reaction; final steady state for a completed pulse)"),
                quantities: first
                    .exits
                    .iter()
                    .enumerate()
                    .filter(|(x, _)| !back.contains(x))
                    .map(|(x, name)| {
                        quantity(name.clone(), temperatures, pressures, &index, |e: &TimeEvolution, _| {
                            last(e).map(|p| 100.0 * p.exit_yields[x] / net(&p))
                        })
                    })
                    .collect(),
            });
        }
    }
    if let Some(r) = reactant {
        // Overall bimolecular-to-bimolecular rates and yields of the reactant (products summed over their
        // channels, the sinks; the return to the reactant excluded).
        let mut targets: Vec<(String, Vec<usize>)> = Vec::new();
        for (x, name) in first.exits.iter().enumerate() {
            let target = match name.split_once("->") {
                Some((_, p)) if p == r => continue,
                Some((_, p)) => format!("{r}->{p}"),
                None if name.starts_with("escape(") => format!("{r}->{name}"),
                None => continue,
            };
            match targets.iter_mut().find(|(n, _)| *n == target) {
                Some((_, members)) => members.push(x),
                None => targets.push((target, vec![x])),
            }
        }
        if !targets.is_empty() {
            let total = |p: &TimePoint, members: &[usize]| {
                members.iter().map(|&x| p.exit_yields[x]).sum::<f64>()
            };
            let net = |p: &TimePoint| targets.iter().map(|(_, m)| total(p, m)).sum::<f64>();
            if let Some(k_inf) = k_inf {
                groups.push(QuantityGroup {
                    title: format!(
                        "Bimolecular-to-bimolecular rate coefficients, overall (chemical activation + thermal), at the last output time: k({r} -> P) = k_inf Y_P (cm^3/s)"
                    ),
                    quantities: targets
                        .iter()
                        .map(|(n, m)| {
                            quantity(n.clone(), temperatures, pressures, &index, |e: &TimeEvolution, t| last(e).map(|p| k_inf(t) * total(&p, m)))
                        })
                        .collect(),
                });
            }
            groups.push(QuantityGroup {
                title: format!("Bimolecular-to-bimolecular yields, overall (chemical activation + thermal), at the last output time (% of the net reaction of {r})"),
                quantities: targets
                    .iter()
                    .map(|(n, m)| {
                        quantity(n.clone(), temperatures, pressures, &index, |e: &TimeEvolution, _| last(e).map(|p| 100.0 * total(&p, m) / net(&p)))
                    })
                    .collect(),
            });
        }
    }
    groups.push(QuantityGroup {
        title: "Total: populations + yields at the last output time (%; 100 for a pulse)".into(),
        quantities: vec![quantity(
            "total",
            temperatures,
            pressures,
            &index,
            |e: &TimeEvolution, _| {
                last(e).map(|p| {
                    100.0
                        * (p.well_populations.iter().sum::<f64>()
                            + p.exit_yields.iter().sum::<f64>())
                })
            },
        )],
    });
    groups
}

/// All groups in the three views: by temperature, by pressure, temperature-pressure tables.
pub fn write_groups<W: Write>(
    out: &mut W,
    temperatures: &[f64],
    pressures: &[f64],
    groups: &[QuantityGroup],
) -> std::io::Result<()> {
    writeln!(
        out,
        "--- Tables by temperature (rows: pressure in torr) ---\n"
    )?;
    for g in groups {
        write_tables_by_temperature(out, temperatures, pressures, g)?;
    }
    writeln!(out, "--- Tables by pressure (rows: temperature in K) ---\n")?;
    for g in groups {
        write_tables_by_pressure(out, temperatures, pressures, g)?;
    }
    writeln!(
        out,
        "--- Temperature-pressure tables (rows: pressure in torr; columns: temperature in K) ---\n"
    )?;
    for g in groups {
        write_temperature_pressure_tables(out, temperatures, pressures, g)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_driver::{
        run_chemical_activation, run_phenomenological_rates, run_thermal_rate_coefficients,
        run_thermal_well_fates, ChemicalActivationRun, SourceSpecification,
    };
    use crate::masterequation::chemical_activation_eigen::{
        EigenSolver, DEFAULT_SUM_RULE_TOLERANCE,
    };
    use crate::masterequation::chemical_activation_network::{
        AbsorbingBarrier, Channel, ChannelDestination, ChemicalActivationOptions, CollisionModel,
        SteadyState,
    };
    use crate::masterequation::chemical_activation_operator::tests::two_well_network;
    use crate::masterequation::chemical_activation_steady_state::LinearSolver;

    const MODEL: CollisionModel = CollisionModel::ExponentialDown {
        cutoff_in_mean_down: 10.0,
    };
    const TEMPERATURES: [f64; 2] = [250.0, 300.0];
    const PRESSURES: [f64; 2] = [10.0, 760.0];

    /// Two wells A, B (B with a sink), entrance A <- R opening at grain 320.
    fn network() -> ChemicalActivationNetwork {
        let mut network = two_well_network();
        network.wells[0].channels.push(network_entrance_channel());
        network
    }

    /// The entrance channel A -> R of `network`.
    fn network_entrance_channel() -> Channel {
        Channel {
            name: "A->reactants".into(),
            destination: ChannelDestination::Products { name: "R".into() },
            threshold_grain: None,
            rate_constant_s_inv: (0..400)
                .map(|i| {
                    if i >= 320 {
                        2.0e7 * ((i - 320) as f64 + 1.0)
                    } else {
                        0.0
                    }
                })
                .collect(),
        }
    }

    fn steady_states(
        network: &ChemicalActivationNetwork,
        steady_state: SteadyState,
    ) -> Vec<ConditionResult> {
        let run = ChemicalActivationRun {
            temperatures_kelvin: TEMPERATURES.to_vec(),
            pressures_torr: PRESSURES.to_vec(),
            options: ChemicalActivationOptions {
                collision_model: MODEL,
                steady_state,
            },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance {
                channels: vec![(0, 2)],
            },
            tolerance: 1e-8,
        };
        run_chemical_activation(network, &run).unwrap()
    }

    fn intermediate() -> SteadyState {
        SteadyState::Intermediate {
            barrier: AbsorbingBarrier::default(),
        }
    }

    fn group<'a>(groups: &'a [QuantityGroup], title_start: &str) -> &'a QuantityGroup {
        groups
            .iter()
            .find(|g| g.title.starts_with(title_start))
            .unwrap_or_else(|| {
                panic!(
                    "no group '{title_start}' in {:?}",
                    groups.iter().map(|g| &g.title).collect::<Vec<_>>()
                )
            })
    }

    /// Sum over the quantities of a group at every condition.
    fn sums(group: &QuantityGroup) -> Vec<f64> {
        let mut out = Vec::new();
        for t in 0..TEMPERATURES.len() {
            for p in 0..PRESSURES.len() {
                out.push(
                    group
                        .quantities
                        .iter()
                        .map(|q| q.values[t][p].unwrap())
                        .sum(),
                );
            }
        }
        out
    }

    fn names(group: &QuantityGroup) -> Vec<&str> {
        group.quantities.iter().map(|q| q.name.as_str()).collect()
    }

    #[test]
    fn energetics_are_in_kcal_per_mol_relative_to_the_reactant() {
        let deck = crate::masterequation::mess_input::parse_mess_input(
            crate::masterequation::chemical_activation_from_mess_input::tests::DECK,
        )
        .unwrap();
        let mut out = Vec::new();
        write_energetics(&mut out, &deck).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.contains("kcal/mol") && text.contains("relative to R"),
            "{text}"
        );
        let row = |name: &str| -> Vec<String> {
            text.lines()
                .find(|l| l.split_whitespace().next() == Some(name))
                .unwrap_or_else(|| panic!("no row {name} in\n{text}"))
                .split_whitespace()
                .map(String::from)
                .collect()
        };
        // W1: ground energy -30, lowest barrier out of it B12 at -5 (B0, barrierless, at the R asymptote 0).
        assert_eq!(row("W1")[1..3], ["-30.00", "-5.00"]);
        assert_eq!(row("R")[1], "0.00");
        assert_eq!(row("B12")[1..4], ["-5.00", "W1", "W2"]);
        assert_eq!(row("B0")[1], "0.00");
    }

    #[test]
    fn steady_state_yields_are_percent_of_the_formed_adducts_and_sum_to_100() {
        let network = network();
        let results = steady_states(&network, intermediate());
        let k_inf = |_t: f64| 2.0e-11;
        let groups = steady_state_groups(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            &results,
            Some("R"),
            Some(&k_inf),
        );
        let yields = group(&groups, "Yields (% of the formed adducts)");
        assert!(
            names(yields).contains(&"A->R")
                && names(yields).contains(&"stab(A)")
                && names(yields).contains(&"escape(B)")
        );
        for s in sums(yields) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        let net = group(&groups, "Yields without the return to R");
        assert!(!names(net).contains(&"A->R"));
        for s in sums(net) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        let k_ca = group(&groups, "Chemical-activation rate coefficients");
        assert!(names(k_ca).contains(&"A->B") && names(k_ca).contains(&"B->A"));
        // Reactant R: bimolecular-to-bimolecular (chemical activation: product P and the escape sink of B) and
        // bimolecular-to-well (stabilization), rates k_inf Phi and yields; capture, return and net reaction.
        let to_bimolecular = group(&groups, "Bimolecular-to-bimolecular rate coefficients");
        assert_eq!(names(to_bimolecular), ["R->P", "R->escape(B)"]);
        let to_wells = group(&groups, "Bimolecular-to-well rate coefficients");
        assert_eq!(names(to_wells), ["R->A", "R->B"]);
        let capture = group(&groups, "Capture, return and net reaction of R");
        assert_eq!(names(capture), ["capture", "return", "net"]);
        let yields_bb = group(&groups, "Bimolecular-to-bimolecular yields");
        let yields_bw = group(&groups, "Bimolecular-to-well yields");
        for ((a, b), t) in sums(yields_bb).iter().zip(sums(yields_bw)).zip(sums(net)) {
            // % of the net reaction: chemically activated products + stabilization = everything not returned.
            assert!(
                (a + b - 100.0).abs() < 1e-6 && (t - 100.0).abs() < 1e-6,
                "{a} + {b}"
            );
        }
        for t in 0..TEMPERATURES.len() {
            for p in 0..PRESSURES.len() {
                let k = 2.0e-11;
                let rate = to_wells.quantities[0].values[t][p].unwrap();
                let phi = group(&groups, "Yields (% of the formed adducts)")
                    .quantities
                    .iter()
                    .find(|q| q.name == "stab(A)")
                    .unwrap()
                    .values[t][p]
                    .unwrap();
                assert!(
                    (rate - k * phi / 100.0).abs() < 1e-12 * k,
                    "k(R->A) = k_inf Phi_stab(A)"
                );
                let net_rate = capture.quantities[2].values[t][p].unwrap();
                let bb = to_bimolecular.quantities[0].values[t][p].unwrap();
                assert!(
                    (100.0 * bb / net_rate - yields_bb.quantities[0].values[t][p].unwrap()).abs()
                        < 1e-8
                );
            }
        }
        // The final steady state gives the overall bimolecular-to-bimolecular rates (no stabilization).
        let total = steady_states(&network, SteadyState::Final);
        let final_groups = steady_state_groups(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            &total,
            Some("R"),
            Some(&k_inf),
        );
        assert_eq!(
            names(group(
                &final_groups,
                "Bimolecular-to-bimolecular rate coefficients, overall"
            )),
            ["R->P", "R->escape(B)"]
        );
        assert!(final_groups
            .iter()
            .all(|g| !g.title.starts_with("Bimolecular-to-well")));
    }

    #[test]
    fn missing_conditions_are_none() {
        let network = network();
        let results = steady_states(&network, intermediate());
        let groups = steady_state_groups(
            &network,
            &[250.0, 300.0, 350.0],
            &PRESSURES,
            &results,
            None,
            None,
        );
        let yields = group(&groups, "Yields (% of the formed adducts)");
        assert!(yields
            .quantities
            .iter()
            .all(|q| q.values[2].iter().all(|v| v.is_none())));
        assert!(yields
            .quantities
            .iter()
            .all(|q| q.values[1].iter().all(|v| v.is_some())));
    }

    #[test]
    fn thermal_branching_and_well_fates_sum_to_100() {
        let network = network();
        let thermal = run_thermal_rate_coefficients(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            MODEL,
            EigenSolver::FullDecomposition,
            DEFAULT_SUM_RULE_TOLERANCE,
        )
        .unwrap();
        let groups = thermal_groups(&network, &TEMPERATURES, &PRESSURES, &thermal);
        assert!(
            names(group(&groups, "Thermal rate coefficients")).starts_with(&["k_uni", "lambda_1"])
        );
        for s in sums(group(&groups, "Thermal branching")) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        group(&groups, "High-pressure rate coefficients");
        group(&groups, "Diagnostics of the thermal eigenpair");
        let fates = run_thermal_well_fates(&network, &TEMPERATURES, &PRESSURES, MODEL).unwrap();
        let fate_groups = well_fate_groups(&network, &TEMPERATURES, &PRESSURES, &fates);
        assert_eq!(fate_groups.len(), 2);
        for g in &fate_groups {
            for s in sums(g) {
                assert!((s - 100.0).abs() < 1e-6, "{}: {s}", g.title);
            }
        }
    }

    #[test]
    fn chemical_activation_and_thermal_contributions_add_up() {
        let network = network();
        let prompt = steady_states(&network, intermediate());
        let total = steady_states(&network, SteadyState::Final);
        let fates = run_thermal_well_fates(&network, &TEMPERATURES, &PRESSURES, MODEL).unwrap();
        let groups = chemical_activation_and_thermal_groups(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            &prompt,
            &fates,
            &total,
        );
        let p = group(&groups, "Chemically activated (prompt)");
        let v = group(&groups, "Through the stabilized wells");
        let s = group(&groups, "Prompt + through the stabilized wells");
        let f = group(&groups, "All together");
        let d = group(&groups, "Difference");
        for (q, quantity) in s.quantities.iter().enumerate() {
            for t in 0..TEMPERATURES.len() {
                for k in 0..PRESSURES.len() {
                    let sum = p.quantities[q].values[t][k].unwrap()
                        + v.quantities[q].values[t][k].unwrap();
                    assert!((quantity.values[t][k].unwrap() - sum).abs() < 1e-12);
                    let diff =
                        quantity.values[t][k].unwrap() - f.quantities[q].values[t][k].unwrap();
                    assert!((d.quantities[q].values[t][k].unwrap() - diff).abs() < 1e-12);
                }
            }
        }
        for total in sums(f) {
            assert!((total - 100.0).abs() < 1e-6, "{total}");
        }
    }

    #[test]
    fn cse_groups_give_the_rates_and_the_yields_of_the_reactant() {
        let network = network();
        let capture = |_t: f64| 2.0e-11;
        let results = run_phenomenological_rates(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            MODEL,
            EigenSolver::FullDecomposition,
            Some(("R", &capture)),
            &crate::masterequation::chemically_significant_eigenvalues::CseMerging::default(),
        )
        .unwrap();
        let groups = cse_groups(&TEMPERATURES, &PRESSURES, &results);
        assert!(names(group(&groups, "Rate coefficients from A")).contains(&"A->B"));
        assert_eq!(
            names(group(
                &groups,
                "Bimolecular-to-bimolecular rate coefficients"
            )),
            ["R->P", "R->escape(B)"]
        );
        assert_eq!(
            names(group(&groups, "Bimolecular-to-well rate coefficients")),
            ["R->A", "R->B"]
        );
        assert_eq!(
            names(group(&groups, "Capture, return and net reaction of R")),
            ["capture", "return", "net"]
        );
        let yields_bb = group(&groups, "Bimolecular-to-bimolecular yields");
        let yields_bw = group(&groups, "Bimolecular-to-well yields");
        for (a, b) in sums(yields_bb).iter().zip(sums(yields_bw)) {
            assert!((a + b - 100.0).abs() < 1e-6, "{a} + {b}");
        }
        for s in sums(group(&groups, "Long-time yields, total")) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        for s in sums(group(&groups, "Thermal fate of A")) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        group(&groups, "Long-time yields, direct");
        group(&groups, "Long-time yields, through the wells");
        let mut out = Vec::new();
        write_cse_species_tables(&mut out, &results).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text.lines().filter(|l| l.starts_with("From\\To")).count(),
            4,
            "{text}"
        );
        assert!(
            text.contains("Temperature = 300 K    Pressure = 760 torr"),
            "{text}"
        );
        assert!(text.contains("chemical eigenvalues"), "{text}");
    }

    #[test]
    fn cse_groups_follow_the_species_by_name_where_wells_are_merged_at_some_conditions() {
        // A and B in fast equilibrium (low isomerization barrier), with the entrance A <- R. At 300 K, 760 Torr a
        // ChemicalEigenvalueMax between Lambda_1 and Lambda_2 merges them into the species A+B (Georgievskii et
        // al. 2013, Sec. IV); at the other conditions they are distinct. Every quantity is looked up by species
        // name and is missing where its species does not exist.
        use crate::masterequation::chemical_activation_operator::assemble_operator;
        use crate::masterequation::chemically_significant_eigenvalues::tests::fast_equilibrium_network;
        use crate::masterequation::chemically_significant_eigenvalues::{phenomenological_rate_coefficients, ChemicalSubspaceCriterion, CseMerging};
        let mut network = fast_equilibrium_network();
        network.wells[0].channels.push(network_entrance_channel());
        let k_capture = 2.0e-11;
        let capture = |_t: f64| k_capture;
        let mut results = run_phenomenological_rates(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            MODEL,
            EigenSolver::FullDecomposition,
            Some(("R", &capture)),
            &CseMerging { chemical_eigenvalue_max: 0.999, ..CseMerging::default() },
        )
        .unwrap();
        let last = results.last_mut().unwrap();
        assert_eq!((last.conditions.temperature_kelvin, last.conditions.pressure_torr), (300.0, 760.0));
        let (l1, l2, l3) = (
            last.rates.chemical_eigenvalues_s_inv[0],
            last.rates.chemical_eigenvalues_s_inv[1],
            last.rates.relaxation_eigenvalue_s_inv,
        );
        assert!(l2 > 10.0 * l1, "the test needs Lambda_2 well above Lambda_1: {l1:e} {l2:e} {l3:e}");
        let merging = CseMerging { chemical_eigenvalue_max: (l1 * l2).sqrt() / l3, criterion: ChemicalSubspaceCriterion::EigenvalueRatio, ..CseMerging::default() };
        let options = ChemicalActivationOptions { collision_model: MODEL, steady_state: SteadyState::Final };
        let op = assemble_operator(&network, &last.conditions, &options).unwrap();
        last.rates =
            phenomenological_rate_coefficients(&network, &op, Some(("R", k_capture)), EigenSolver::FullDecomposition, &merging)
                .unwrap();
        assert_eq!(last.rates.wells, ["A+B"]);

        let groups = cse_groups(&TEMPERATURES, &PRESSURES, &results);
        let values = |title: &str, name: &str| -> Vec<Option<f64>> {
            let g = group(&groups, title);
            let q = g.quantities.iter().find(|q| q.name == name).unwrap_or_else(|| panic!("no '{name}' in {:?}", names(g)));
            q.values.iter().flatten().copied().collect()
        };
        let present = |v: Vec<Option<f64>>| v.iter().map(|x| x.is_some()).collect::<Vec<_>>();
        let only_last = [false, false, false, true];
        let all_but_last = [true, true, true, false];
        assert_eq!(present(values("Rate coefficients from A ", "A->B")), all_but_last);
        assert_eq!(present(values("Rate coefficients from A+B", "A+B loss")), only_last);
        assert_eq!(present(values("Rate coefficients from A+B", "A+B->P")), only_last);
        assert_eq!(names(group(&groups, "Bimolecular-to-well rate coefficients")), ["R->A", "R->B", "R->A+B"]);
        assert_eq!(present(values("Bimolecular-to-well rate coefficients", "R->A")), all_but_last);
        assert_eq!(present(values("Bimolecular-to-well rate coefficients", "R->A+B")), only_last);
        // The merged value is the one of the merged species.
        let k_r_ab = values("Bimolecular-to-well rate coefficients", "R->A+B")[3].unwrap();
        assert!((k_r_ab / last_rates(&results).reactant.as_ref().unwrap().to_well_cm3_s[0] - 1.0).abs() < 1e-12);
        // The prompt branching and the long-time yields close at every condition.
        let sum_present = |title: &str| -> Vec<f64> {
            let g = group(&groups, title);
            (0..TEMPERATURES.len() * PRESSURES.len())
                .map(|c| g.quantities.iter().filter_map(|q| q.values[c / PRESSURES.len()][c % PRESSURES.len()]).sum())
                .collect()
        };
        for (a, b) in sum_present("Bimolecular-to-bimolecular yields").iter().zip(sum_present("Bimolecular-to-well yields")) {
            assert!((a + b - 100.0).abs() < 1e-6, "{a} + {b}");
        }
        for s in sum_present("Long-time yields, total") {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        assert!((sum_present("Thermal fate of A+B")[3] - 100.0).abs() < 1e-6);
        // The number of species at every condition.
        let species: Vec<f64> = values("Diagnostics of the CSE solution", "species").into_iter().map(Option::unwrap).collect();
        assert_eq!(species, [2.0, 2.0, 2.0, 1.0]);
    }

    fn last_rates(
        results: &[PhenomenologicalConditionResult],
    ) -> &crate::masterequation::chemically_significant_eigenvalues::PhenomenologicalRates {
        &results.last().unwrap().rates
    }

    #[test]
    fn the_network_summary_lists_wells_channels_models_and_sinks() {
        use crate::masterequation::chemical_activation_from_mess_input::{
            chemical_activation_model_from_mess, MessNetworkSettings,
        };
        let deck = crate::masterequation::mess_input::parse_mess_input(
            crate::masterequation::chemical_activation_from_mess_input::tests::DECK,
        )
        .unwrap();
        let model =
            chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap();
        let mut out = Vec::new();
        write_network_summary(&mut out, &deck, &model.network, &model.entrance_channels).unwrap();
        let text = String::from_utf8(out).unwrap();
        let row = |name: &str| {
            text.lines()
                .find(|l| l.split_whitespace().next() == Some(name))
                .unwrap_or_else(|| panic!("{name}\n{text}"))
        };
        assert!(row("W1").contains("<dE_down>"), "{text}");
        assert!(
            row("B12").contains("W1")
                && row("B12").contains("W2")
                && row("B12").contains("rigid transition state"),
            "{text}"
        );
        assert!(
            row("B0").contains("inverse Laplace transform") && row("B0").contains("entrance"),
            "{text}"
        );
        assert!(row("B2P").contains("P"), "{text}");
    }

    #[test]
    fn time_integration_tables_and_groups() {
        use crate::masterequation::chemical_activation_sources::thermal_entrance_source;
        use crate::masterequation::direct_time_integration::{
            integrate_master_equation, InitialState, TimeIntegrationSettings,
        };
        let network = network();
        let mut evolutions = Vec::new();
        for &t in &TEMPERATURES {
            for &p in &PRESSURES {
                let f = thermal_entrance_source(&network, &[(0, 2)], crate::constants::KB_CM * t)
                    .unwrap();
                let settings = TimeIntegrationSettings {
                    times_s: vec![1e-9, 1e-6, 1e-3, 1.0],
                    ..Default::default()
                };
                let options = ChemicalActivationOptions {
                    collision_model: MODEL,
                    steady_state: SteadyState::Final,
                };
                let conditions = crate::masterequation::chemical_activation_network::Conditions {
                    temperature_kelvin: t,
                    pressure_torr: p,
                };
                evolutions.push(
                    integrate_master_equation(
                        &network,
                        &conditions,
                        &options,
                        &f,
                        InitialState::Pulse,
                        &settings,
                    )
                    .unwrap(),
                );
            }
        }
        let mut out = Vec::new();
        write_time_evolution_tables(&mut out, &network, &evolutions).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text.lines()
                .filter(|l| l.trim_start().starts_with("t(s)"))
                .count(),
            4,
            "{text}"
        );
        assert!(
            text.contains("Temperature = 300 K    Pressure = 760 torr"),
            "{text}"
        );
        let k_inf = |_t: f64| 2.0e-11;
        let groups = time_integration_groups(
            &network,
            &TEMPERATURES,
            &PRESSURES,
            &evolutions,
            Some("R"),
            Some(&k_inf),
        );
        // Reactant: overall bimolecular-to-bimolecular rates and yields (chemical activation + thermal) at the
        // last output time, products summed over their channels.
        assert_eq!(
            names(group(
                &groups,
                "Bimolecular-to-bimolecular rate coefficients, overall"
            )),
            ["R->P", "R->escape(B)"]
        );
        for s in sums(group(&groups, "Bimolecular-to-bimolecular yields, overall")) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        for s in sums(group(&groups, "Yields without the return to R")) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
        let yields = group(&groups, "Yields at the last output time");
        let left = group(&groups, "Populations left in the wells");
        for (a, b) in sums(yields).iter().zip(sums(left)) {
            assert!((a + b - 100.0).abs() < 1e-6, "{a} + {b}");
        }
        for s in sums(group(&groups, "Total")) {
            assert!((s - 100.0).abs() < 1e-6, "{s}");
        }
    }

    #[test]
    fn groups_are_written_in_three_views() {
        let network = network();
        let results = steady_states(&network, intermediate());
        let groups = steady_state_groups(&network, &TEMPERATURES, &PRESSURES, &results, None, None);
        let mut out = Vec::new();
        write_groups(&mut out, &TEMPERATURES, &PRESSURES, &groups).unwrap();
        let text = String::from_utf8(out).unwrap();
        for view in [
            "by temperature:",
            "by pressure:",
            "temperature-pressure tables:",
        ] {
            assert!(
                text.contains(&format!("Yields (% of the formed adducts), {view}")),
                "{view}"
            );
        }
    }
}
