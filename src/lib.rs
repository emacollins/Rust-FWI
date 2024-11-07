pub mod utils;

pub fn ffmc(temp: f64, mut rh: f64, wind: f64, precip: f64, ffmc_prev: f64) -> f64 {

    let mut mo: f64; // Fine fuel moisture content from previous day
    let ed: f64;
    let m: f64;

    rh = utils::min(100.0, rh);

    mo = 147.2 * (101.0 - ffmc_prev) / (59.5 + ffmc_prev);

    if precip > 0.5 {
        let mut mr: f64;
        let rf: f64;

        rf = precip - 0.5;

        if mo <= 150.0 {
            mr = mo + 42.5 * rf * (-100.0 / (251.1 - mo)).exp() * (1.0-(-6.93 / rf).exp());
        } else {
            mr = mo + 42.5 * rf * (-100.0 / (251.1 - mo)).exp() * (1.0-(-6.93 / rf).exp()) + (0.0015 * (mo - 150.0).powi(2) * rf.sqrt());
        }

        if mr > 250.0 {
            mr = 250.0;
        }
        mo = mr;
    }

    ed = 0.942 * rh.powf(0.679) + 11.0 * ((rh - 100.0) / 10.0).exp() + 0.18 * (21.1 - temp) * (1.0 - (-0.115 * rh).exp());

    if mo > ed {
        let ko: f64;
        let kd: f64;

        ko = 0.424 * (1.0 - (rh / 100.0).powf(1.7)) + 0.0694 * wind.sqrt() * (1.0 - (rh / 100.0).powi(8));
        kd = ko * 0.581 * (0.0365 * temp).exp();
        m = ed + (mo - ed) * 10.0_f64.powf(-kd);

    } else {
        let ew: f64;

        ew = 0.618 * rh.powf(0.753) + 10.0 * ((rh - 100.0) / 10.0) + 0.18 * (21.1 - temp) * (1.0 - (-0.115 * rh).exp());

        if mo < ew {

            let k1: f64;
            let kw: f64;

            k1 = 0.424 * (1.0 - (100.0 - rh / 100.0).powf(1.7)) + 0.0694 * wind.sqrt() * (1.0 - ((100.0 - rh) / 100.0).powi(8));

            kw = k1 * 0.581 * (0.0365 * temp).exp();

            m = ew - (ew - mo) * 10.0_f64.powf(-kw);
        } else {
            m = mo;
        }
    }

    59.5 * (250.0 - m) / (147.2 + m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ffmc() {
        let ffmc: f64 = ffmc(17.0, 42.0, 25.0, 0.0, 85.0);
        assert_eq!(ffmc, 87.692980092774448);
    }
}