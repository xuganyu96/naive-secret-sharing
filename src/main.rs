use std::fmt::write;

use secretsharing::f2x::*;

fn main() {
    for limb in 1..=0xFFFF {
        let lhs = F2x::<8>::from_limbs([0, 0, 0, 0, 0, 0, 0, limb]);
        let rhs = WideF2x::from_f2x(
            F2x::<8>::ONE,
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x2000, 0x0000, 0x0008, 0x0000, 0x0801,
            ]),
        );
        let gcd = WideF2x::gcd(&lhs.widen(), &rhs);
        if gcd != WideF2x::<8>::ONE {
            println!("{limb:X}");
        }
    }
}
