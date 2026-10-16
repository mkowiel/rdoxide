//! Stereo descriptors (upstream `CIPLabeler/Descriptor.h`).

/// `CIPLabeler::Descriptor`. The lower case descriptors of upstream (`r`,
/// `s`, `m`, `p`) are the `Pseudo*` variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Descriptor {
    /// unspecified
    None,
    Unknown,
    /// other
    Ns,
    R,
    S,
    PseudoR,
    PseudoS,
    SeqTrans,
    SeqCis,
    E,
    Z,
    M,
    P,
    PseudoM,
    PseudoP,
    Sp4,
    Tbpy5,
    Oc6,
}

impl Descriptor {
    /// `to_string(Descriptor)`
    pub fn as_str(self) -> &'static str {
        match self {
            Descriptor::None => "NONE",
            Descriptor::Unknown => "UNKNOWN",
            Descriptor::Ns => "ns",
            Descriptor::R => "R",
            Descriptor::S => "S",
            Descriptor::PseudoR => "r",
            Descriptor::PseudoS => "s",
            Descriptor::SeqTrans => "e",
            Descriptor::SeqCis => "z",
            Descriptor::E => "E",
            Descriptor::Z => "Z",
            Descriptor::M => "M",
            Descriptor::P => "P",
            Descriptor::PseudoM => "m",
            Descriptor::PseudoP => "p",
            Descriptor::Sp4 => "SP_4",
            Descriptor::Tbpy5 => "TBPY_5",
            Descriptor::Oc6 => "OC_6",
        }
    }
}
