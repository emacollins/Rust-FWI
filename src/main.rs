extern crate rust_fwi as fwi;


fn main() {
   let ffmc: f64 = fwi::ffmc(17.0, 42.0, 25.0, 0.0, 85.0);
   println!("ffmc = {ffmc}");

}
