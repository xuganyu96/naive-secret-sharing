use secretsharing::*;

fn main() {
    let (high, low) = widening_clmul(0xE223, 0x672F);
    println!("{high:04X}, {low:04X}");
}
