use super::reaction_network::WellDefinition;

/// Maps (well_index, local_grain_index) to global stacked index.
pub struct GlobalLayout {
    pub total_state_count: usize,
    start_index_per_well: Vec<usize>,
    lowest_included_per_well: Vec<usize>,
    one_past_highest_included_per_well: Vec<usize>,
}

impl GlobalLayout {
    pub(crate) fn from_wells(wells: &[WellDefinition]) -> Result<Self, String> {
        let mut start_index_per_well = Vec::with_capacity(wells.len());
        let mut lowest = Vec::with_capacity(wells.len());
        let mut highest = Vec::with_capacity(wells.len());

        let mut running = 0usize;
        for w in wells {
            if w.one_past_highest_included_grain_index <= w.lowest_included_grain_index {
                return Err(format!("Invalid grain bounds in well {}", w.well_name));
            }
            start_index_per_well.push(running);
            lowest.push(w.lowest_included_grain_index);
            highest.push(w.one_past_highest_included_grain_index);

            running += w.one_past_highest_included_grain_index - w.lowest_included_grain_index;
        }

        Ok(Self {
            total_state_count: running,
            start_index_per_well,
            lowest_included_per_well: lowest,
            one_past_highest_included_per_well: highest,
        })
    }

    pub(crate) fn global_index_of(
        &self,
        well_index: usize,
        local_grain_index: usize,
    ) -> Result<usize, String> {
        if well_index >= self.start_index_per_well.len() {
            return Err("Well index out of bounds".into());
        }
        if local_grain_index < self.lowest_included_per_well[well_index]
            || local_grain_index >= self.one_past_highest_included_per_well[well_index]
        {
            return Err("Local grain index out of included bounds".into());
        }
        let offset = local_grain_index - self.lowest_included_per_well[well_index];
        Ok(self.start_index_per_well[well_index] + offset)
    }
}
