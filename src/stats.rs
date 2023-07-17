use std::{collections::BTreeMap, fmt::Display, sync::Mutex};

use lazy_static::lazy_static;

lazy_static! {
    pub static ref STATS: Mutex<ProtocolStats> = Mutex::new(ProtocolStats {
        time: BTreeMap::new(),
        total_time: 0.0,
    });
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct ProtocolStats {
    pub time: BTreeMap<ProtocolLowOps, (usize, f32)>,
    pub total_time: f32,
}

impl Display for ProtocolStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Protocol Stats {{")?;
        for (op, (occurance, time)) in &self.time {
            writeln!(
                f,
                "{}: x{}, {:.2}s ({:.2}%)",
                op,
                occurance,
                time,
                time * 100. / self.total_time
            )?;
        }
        writeln!(f, "Total Time: {:.2}s", self.total_time)?;
        writeln!(f, "}}")?;

        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ProtocolLowOps {
    AddMod,
    SubMod,
    DoubleMod,
    MulMod,
    SquareMod,
    InverseMod,
}

impl Display for ProtocolLowOps {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProtocolLowOps::AddMod => f.write_str("Add Mod"),
            ProtocolLowOps::SubMod => f.write_str("Sub Mod"),
            ProtocolLowOps::DoubleMod => f.write_str("Double Mod"),
            ProtocolLowOps::MulMod => f.write_str("Mul Mod"),
            ProtocolLowOps::SquareMod => f.write_str("Square Mod"),
            ProtocolLowOps::InverseMod => f.write_str("Inverse Mod"),
        }
    }
}

impl ProtocolStats {
    pub fn add_time(op: ProtocolLowOps, time: f32) {
        let mut stat = STATS.lock().unwrap();
        let entry = stat.time.entry(op).or_insert((0, 0.0));
        *entry = (entry.0 + 1, entry.1 + time);
        stat.total_time += time;
    }

    pub fn reset() {
        let mut stat = STATS.lock().unwrap();
        stat.time.clear();
    }

    pub fn stats() -> Self {
        STATS.lock().unwrap().clone()
    }
}
