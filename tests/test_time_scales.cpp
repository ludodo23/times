#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_approx.hpp>

#include "time_scales.hpp"

using namespace ts;
using Catch::Matchers::WithinAbs;


static constexpr double EPS = 1e-12;

TEST_CASE("JulianDate normalization", "[jd]")
{
    JulianDate jd {2451545.0, 1.25};

    auto n = jd.normalized();

    REQUIRE(n.jd1 == Catch::Approx(2451546.0));
    REQUIRE(n.jd2 == Catch::Approx(0.25));
}

TEST_CASE("JulianDate subtraction preserves precision", "[jd]")
{
    JulianDate a {2451545.0, 1e-9};
    JulianDate b {2451545.0, 0.0};

    double d = a - b;

    REQUIRE(d == Catch::Approx(1e-9));
}

TEST_CASE("Calendar to JD J2000", "[calendar]")
{
    auto jd = calendar_to_jd(2000, 1, 1, 12, 0, 0.0);

    REQUIRE(jd.total() == Catch::Approx(2451545.0).epsilon(1e-15));
}

TEST_CASE("JD to calendar J2000", "[calendar]")
{
    auto cal = jd_to_calendar(JulianDate{2451545.0, 0.0});

    REQUIRE(cal.year == 2000);
    REQUIRE(cal.month == 1);
    REQUIRE(cal.day == 1);

    REQUIRE(cal.hour == 12);
    REQUIRE(cal.minute == 0);

    REQUIRE(cal.second == Catch::Approx(0.0));
}

TEST_CASE("Modified Julian Date conversion", "[mjd]")
{
    JulianDate jd {2451545.0, 0.0};

    auto mjd = jd_to_mjd(jd);

    REQUIRE(mjd.total() == Catch::Approx(51544.5));
}

TEST_CASE("MJD roundtrip", "[mjd]")
{
    JulianDate jd {2455197.5, 0.123456789};

    auto mjd = jd_to_mjd(jd);
    auto back = mjd_to_jd(mjd);

    REQUIRE(back.total() == Catch::Approx(jd.total()).epsilon(1e-15));
}

TEST_CASE("UTC to TAI", "[utc][tai]")
{
    JulianDate utc = calendar_to_jd(2017, 1, 1, 0, 0, 0);

    auto tai = utc_to_tai(utc);

    double dt = (tai - utc) * SEC_PER_DAY;

    REQUIRE(dt == Catch::Approx(37.0).margin(1e-9));
}

TEST_CASE("TAI to TT", "[tai][tt]")
{
    JulianDate tai {2457754.5, 0.0};

    auto tt = tai_to_tt(tai);

    double dt = (tt - tai) * SEC_PER_DAY;

    REQUIRE(dt == Catch::Approx(32.184).margin(1e-12));
}

TEST_CASE("UTC -> TAI -> UTC roundtrip", "[roundtrip]")
{
    JulianDate utc = calendar_to_jd(2024, 5, 1, 12, 30, 15.0);

    auto tai = utc_to_tai(utc);
    auto back = tai_to_utc(tai);

    REQUIRE(back.total() == Catch::Approx(utc.total()).epsilon(1e-14));
}

TEST_CASE("TT to TCG monotonicity", "[tcg]")
{
    JulianDate tt1 = calendar_to_jd(2024, 1, 1, 0, 0, 0);
    JulianDate tt2 = calendar_to_jd(2034, 1, 1, 0, 0, 0);

    auto tcg1 = tt_to_tcg(tt1);
    auto tcg2 = tt_to_tcg(tt2);

    REQUIRE(tcg2.total() > tcg1.total());
}

TEST_CASE("TDB minus TT stays small", "[tdb]")
{
    JulianDate tt = calendar_to_jd(2025, 1, 1, 0, 0, 0);

    double delta = tdb_minus_tt(tt);

    REQUIRE(std::abs(delta) < 0.002);
}

TEST_CASE("Leap second boundary", "[leapsecond]")
{
    auto before = calendar_to_jd(CalendarDate(2016, 12, 31, 23, 59, 59.0));
    auto after  = calendar_to_jd(CalendarDate(2017, 1, 1, 0, 0, 0.0));

    REQUIRE(after.total() > before.total());
}

TEST_CASE("High precision fractional seconds", "[precision]")
{
    auto jd = calendar_to_jd(CalendarDate(
        2024,
        1,
        1,
        0,
        0,
        0.123456789
    ));

    auto cal = jd_to_calendar(jd);

    REQUIRE(cal.second == Catch::Approx(0.123456789).epsilon(1e-12));
}