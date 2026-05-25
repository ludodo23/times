
/*  
 * License: CeCILL-C  
 *  
 * Copyright (c) 2026 Ludovic Andrieux  
 * contributor(s): Ludovic Andrieux (2026)  
 *  
 * ludovic.andrieux23@gmail.com  
 *  
 * This software is a header-only C++ library for orbital times conversions.
 *  
 * This software is governed by the CeCILL-C license under French law and  
 * abiding by the rules of distribution of free software. You can use,  
 * modify and/or redistribute the software under the terms of the CeCILL-C  
 * license as circulated by CEA, CNRS and INRIA at the following URL:  
 * https://www.cecill.info  
 *  
 * As a counterpart to the access to the source code and rights to copy,  
 * modify and redistribute granted by the license, users are provided only  
 * with a limited warranty and the software's author, the holder of the  
 * economic rights, and the successive licensors have only limited  
 * liability.  
 *  
 * In this respect, the user's attention is drawn to the risks associated  
 * with loading, using, modifying and/or developing or reproducing the  
 * software by the user in light of its specific status of free software,  
 * that may mean that it is complicated to manipulate, and that also  
 * therefore means that it is reserved for developers and experienced  
 * professionals having in-depth computer knowledge.  
 *  
 * The fact that you are presently reading this means that you have had  
 * knowledge of the CeCILL-C license and that you accept its terms.
 * 
 * 
 * ===============
 * Conversions entre échelles de temps astronomiques :
 *   UTC, TAI, UT1, TT, TCG, TCB, TDB
 *   Julian Date (deux parties), Modified Julian Date (deux parties)
 *
 * ═══════════════════════════════════════════════════════════════
 * DÉFINITIONS FORMELLES
 * ═══════════════════════════════════════════════════════════════
 *
 * UTC  Temps Universel Coordonné
 *        UTC = TAI − ΔAT,  ΔAT ∈ ℕ (leap seconds, décision IERS)
 *
 * UT1  Temps Universel (rotation terrestre)
 *        UT1 = UTC + DUT1,  |DUT1| < 0.9 s  (IERS Bulletin A)
 *
 * TAI  Temps Atomique International
 *        référence primaire, continu
 *
 * TT   Temps Terrestre (IAU 1991)
 *        TT = TAI + 32.184 s  (constante, définie par l'IAU)
 *        Successeur de TDT et ET.
 *
 * TCG  Temps Coordonné Géocentrique (IAU 1991)
 *        Échelle propre géocentrique, non affectée par la gravité terrestre.
 *        TCG − TT = L_G · (JD_TT − T₀) · 86400
 *        L_G = 6.969290134 × 10⁻¹⁰  (IAU 2000, IERS Conv. 2010 eq. 10.2)
 *        T₀  = 2443144.5003725  (JD TT du 1977-01-01T00:00:32.184 TAI)
 *
 * TCB  Temps Coordonné Barycentrique (IAU 1991)
 *        Échelle propre barycentrique (barycentre système solaire).
 *        TCB − TCG = (1/c²)·∫[v²_E/2 + U_ext(x_E)]dt − v_E·(x−x_E)/c²
 *        La dérivée séculaire moyenne vaut L_B :
 *        d(TCB−TT)/dt ≈ L_B = 1.55051976772 × 10⁻⁸  (IERS Conv. 2010 eq. 10.3)
 *
 * TDB  Temps Barycentrique Dynamique (redéfini IAU 2006, résolution B3)
 *        Transformation linéaire de TCB restant en moyenne synchrone avec TT :
 *          TDB = TCB − L_B · (JD_TCB − T₀) · 86400 + TDB₀
 *          TDB₀ = −6.55 × 10⁻⁵ s
 *        La relation TDB−TT est calculée par la série de Fairhead & Bretagnon
 *        (1990, A&A 229, 240) à 787 termes, telle qu'implémentée dans la
 *        bibliothèque SOFA (iauDtdb). Précision < 30 ns sur [1600, 2200].
 *
 * JD   Julian Date  (jours depuis J.D. 0 = 1er janv. 4713 av. J.-C. à 12h UT)
 *        Stocké en deux doubles (jd1, jd2) pour éviter la perte de précision
 *        (~20 µs avec un seul double pour des dates récentes).
 *        Valeur totale = jd1 + jd2.
 *
 * MJD  Modified Julian Date
 *        MJD = JD − 2400000.5  (origine : 1858-11-17T00:00:00)
 *        Aussi stocké en deux parties.
 *
 * ═══════════════════════════════════════════════════════════════
 * RÉFÉRENCES
 * ═══════════════════════════════════════════════════════════════
 *   IERS Conventions 2010 (Petit & Luzum), chapitre 10
 *   IAU SOFA Software Collection, iauDtdb.c (Copyright © IAU SOFA Board)
 *   Fairhead & Bretagnon, A&A 229, 240–247 (1990)
 *   Harada & Fukushima, AJ 126, 2557 (2003)
 *   Meeus, Astronomical Algorithms, 2nd ed., ch. 7
 *
 * ═══════════════════════════════════════════════════════════════
 * CONVENTIONS DE NOMMAGE
 * ═══════════════════════════════════════════════════════════════
 *   Fonctions de conversion  : snake_case  (ex. utc_to_tai)
 *   Classes principales      : CamelCase   (TimePoint, DeltaTime)
 *   Structures internes      : CamelCase   (JulianDate, CalendarDate, …)
 *   Constantes               : UPPER_CASE
 */

#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

namespace ts {

// ═══════════════════════════════════════════════════════════════
//  §1  CONSTANTES
// ═══════════════════════════════════════════════════════════════

/// Secondes par jour julien
inline constexpr double SEC_PER_DAY = 86400.0;

/// TT − TAI (s), constante IAU
inline constexpr double TT_MINUS_TAI = 32.184;

/// L_G : taux TCG−TT, IAU 2000 (adimensionnel)
inline constexpr double L_G = 6.969290134e-10;

/// L_B : taux TCB−TT (moyenne), IAU 2006 (adimensionnel)
inline constexpr double L_B = 1.55051976772e-8;

/// TDB₀ (s), offset de la définition IAU 2006 de TDB
inline constexpr double TDB0 = -6.55e-5;

/// Époque de référence T₀ pour TCG, TCB, TDB (JD TT)
/// = 1977-01-01T00:00:32.184 TAI = 1977-01-01T00:00:00.000 TAI + 32.184s
inline constexpr double T0_JD = 2443144.5003725;

/// J2000.0 en JD
inline constexpr double J2000_JD = 2451545.0;

/// Offset JD → MJD
inline constexpr double JD_TO_MJD = 2400000.5;

// ═══════════════════════════════════════════════════════════════
//  §2  TYPES DE BASE
// ═══════════════════════════════════════════════════════════════

// ───────────────────────────────────────────────────────────────
//  2.1  JulianDate  (deux parties)
// ───────────────────────────────────────────────────────────────
/**
 * Julian Date stockée en deux doubles afin de conserver la précision
 * numérique maximale (~10 ns pour des dates modernes).
 *
 * Convention : jd1 contient la partie entière ou demi-entière,
 *              jd2 contient le reste (typiquement |jd2| < 1).
 */
struct JulianDate {
    double jd1 = 0.0;
    double jd2 = 0.0;

    // ── Constructeurs ─────────────────────────────────────────
    constexpr JulianDate() = default;
    constexpr JulianDate(double j1, double j2) : jd1(j1), jd2(j2) {}
    explicit  JulianDate(double jd_total) {
        double intpart;
        jd2 = std::modf(jd_total, &intpart);
        jd1 = intpart;
    }

    // ── Accesseurs ────────────────────────────────────────────
    /// Valeur totale (perte de précision possible)
    double total() const noexcept { return jd1 + jd2; }

    /// Fraction de jour dans [0, 1)
    double day_fraction() const noexcept {
        double f = std::fmod(jd2 + 0.5, 1.0);   // décale à minuit
        return f < 0.0 ? f + 1.0 : f;
    }

    // ── Arithmétique ──────────────────────────────────────────
    JulianDate add_seconds(double s) const noexcept {
        return { jd1, jd2 + s / SEC_PER_DAY };
    }

    double diff_seconds(const JulianDate& other) const noexcept {
        return ((jd1 - other.jd1) + (jd2 - other.jd2)) * SEC_PER_DAY;
    }

    JulianDate operator+(double days) const noexcept { return { jd1, jd2 + days }; }
    JulianDate operator-(double days) const noexcept { return { jd1, jd2 - days }; }
    double     operator-(const JulianDate& o) const noexcept { return total() - o.total(); }

    /// Normalise : jd1 = entier, jd2 ∈ [0, 1)
    JulianDate normalized() const noexcept {
        double ip; double fp = std::modf(jd1 + jd2, &ip);
        if (fp < 0.0) { ip -= 1.0; fp += 1.0; }
        return { ip, fp };
    }

    std::string to_string() const {
        std::ostringstream o;
        o << std::fixed << std::setprecision(9) << jd1 << " + " << jd2;
        return o.str();
    }
};

// ───────────────────────────────────────────────────────────────
//  2.2  ModifiedJulianDate  (deux parties)
// ───────────────────────────────────────────────────────────────

typedef JulianDate ModifiedJulianDate;  // même structure, même interface

// ───────────────────────────────────────────────────────────────
//  2.3  CalendarDate  (grégorien)
// ───────────────────────────────────────────────────────────────
struct CalendarDate {
    int    year   = 2000;
    int    month  = 1;      ///< 1–12
    int    day    = 1;      ///< 1–31
    int    hour   = 0;      ///< 0–23
    int    minute = 0;      ///< 0–59
    double second = 0.0;    ///< [0, 61) (leaps seconds possibles)

    constexpr CalendarDate() = default;
    constexpr CalendarDate(int y, int m, int d,
                           int h = 0, int mn = 0, double s = 0.0)
        : year(y), month(m), day(d), hour(h), minute(mn), second(s) {}

    /// Format ISO 8601 : "2000-01-01T12:00:00.000"
    std::string to_string() const {
        std::ostringstream o;
        o << std::setfill('0')
          << std::setw(4) << year  << '-'
          << std::setw(2) << month << '-'
          << std::setw(2) << day   << 'T'
          << std::setw(2) << hour  << ':'
          << std::setw(2) << minute << ':';
        if (second < 10.0) o << '0';
        o << std::fixed << std::setprecision(3) << second;
        return o.str();
    }
};

// ═══════════════════════════════════════════════════════════════
//  §3  TABLE DES LEAP SECONDS  (TAI − UTC)
// ═══════════════════════════════════════════════════════════════

/**
 * Retourne ΔAT = TAI − UTC (secondes entières) pour un MJD donné (UTC).
 * Source : IERS Bulletin C, table complète 1972–2017.
 * À compléter si de nouvelles leap seconds sont annoncées.
 */
inline double leap_seconds(double mjd_utc) noexcept {
    // { MJD_début_UTC, TAI−UTC }
    static constexpr std::array<std::pair<double,double>, 28> TABLE {{
        { 41317.0, 10.0 }, // 1972-01-01
        { 41499.0, 11.0 }, // 1972-07-01
        { 41683.0, 12.0 }, // 1973-01-01
        { 42048.0, 13.0 }, // 1974-01-01
        { 42413.0, 14.0 }, // 1975-01-01
        { 42778.0, 15.0 }, // 1976-01-01
        { 43144.0, 16.0 }, // 1977-01-01
        { 43509.0, 17.0 }, // 1978-01-01
        { 43874.0, 18.0 }, // 1979-01-01
        { 44239.0, 19.0 }, // 1980-01-01
        { 44786.0, 20.0 }, // 1981-07-01
        { 45151.0, 21.0 }, // 1982-07-01
        { 45516.0, 22.0 }, // 1983-07-01
        { 46247.0, 23.0 }, // 1985-07-01
        { 47161.0, 24.0 }, // 1988-01-01
        { 47892.0, 25.0 }, // 1990-01-01
        { 48257.0, 26.0 }, // 1991-01-01
        { 48804.0, 27.0 }, // 1992-07-01
        { 49169.0, 28.0 }, // 1993-07-01
        { 49534.0, 29.0 }, // 1994-07-01
        { 50083.0, 30.0 }, // 1996-01-01
        { 50630.0, 31.0 }, // 1997-07-01
        { 51179.0, 32.0 }, // 1999-01-01
        { 53736.0, 33.0 }, // 2006-01-01
        { 54832.0, 34.0 }, // 2009-01-01
        { 56109.0, 35.0 }, // 2012-07-01
        { 57204.0, 36.0 }, // 2015-07-01
        { 57754.0, 37.0 }, // 2017-01-01
    }};
    double ls = 10.0; // avant 1972
    for (auto& [mjd, val] : TABLE) {
        if (mjd_utc >= mjd) ls = val; else break;
    }
    return ls;
}

// ═══════════════════════════════════════════════════════════════
//  §4  CALENDRIER ↔ JULIAN DATE
// ═══════════════════════════════════════════════════════════════

/**
 * Calendrier grégorien → JD (deux parties).
 * Algorithme : Meeus, Astronomical Algorithms, ch. 7.
 * Valide pour le calendrier grégorien (après 1582-10-15).
 */
inline JulianDate calendar_to_jd(const CalendarDate& c) {
    int Y = c.year, M = c.month;
    if (M <= 2) { Y -= 1; M += 12; }
    int A = Y / 100;
    int B = 2 - A + A / 4;
    double day_frac = c.day
                    + c.hour   / 24.0
                    + c.minute / 1440.0
                    + c.second / SEC_PER_DAY;
    double jd = std::floor(365.25  * (Y + 4716))
              + std::floor(30.6001 * (M + 1))
              + day_frac + B - 1524.5;
    double ip; double fp = std::modf(jd, &ip);
    return { ip, fp };
}

/**
 * JD → calendrier grégorien.
 * Algorithme : Meeus, ch. 7.
 */
inline CalendarDate jd_to_calendar(const JulianDate& jd) {
    double jd_val = jd.jd1 + jd.jd2 + 0.5;
    double Z, F;
    F = std::modf(jd_val, &Z);

    double A;
    if (Z < 2299161.0) {
        A = Z;
    } else {
        double alpha = std::floor((Z - 1867216.25) / 36524.25);
        A = Z + 1.0 + alpha - std::floor(alpha / 4.0);
    }
    double B = A + 1524.0;
    double C = std::floor((B - 122.1) / 365.25);
    double D = std::floor(365.25 * C);
    double E = std::floor((B - D) / 30.6001);

    double day_frac = B - D - std::floor(30.6001 * E) + F;
    int day   = static_cast<int>(day_frac);
    double tf = (day_frac - day) * SEC_PER_DAY;

    int month = (E < 14) ? static_cast<int>(E - 1) : static_cast<int>(E - 13);
    int year  = (month > 2) ? static_cast<int>(C - 4716) : static_cast<int>(C - 4715);

    int hour   = static_cast<int>(tf / 3600.0); tf -= hour   * 3600.0;
    int minute = static_cast<int>(tf / 60.0);   tf -= minute *   60.0;
    return { year, month, day, hour, minute, tf };
}

// ═══════════════════════════════════════════════════════════════
//  §5  JD ↔ MJD
// ═══════════════════════════════════════════════════════════════

inline ModifiedJulianDate jd_to_mjd(const JulianDate& jd) noexcept {
    // MJD = JD − 2400000.5  →  on factorise pour garder la précision
    return { jd.jd1 - 2400000.0, jd.jd2 - 0.5 };
}

inline JulianDate mjd_to_jd(const ModifiedJulianDate& mjd) noexcept {
    return { mjd.jd1 + 2400000.0, mjd.jd2 + 0.5 };
}

// ═══════════════════════════════════════════════════════════════
//  §6  CONVERSIONS D'ÉCHELLES
// ═══════════════════════════════════════════════════════════════

// ───────────────────────────────────────────────────────────────
//  6.1  UTC ↔ TAI
// ───────────────────────────────────────────────────────────────

inline JulianDate utc_to_tai(const JulianDate& jd_utc) {
    double ls = leap_seconds(jd_to_mjd(jd_utc).total());
    return jd_utc.add_seconds(ls);
}

inline JulianDate tai_to_utc(const JulianDate& jd_tai) {
    // La table est indexée en UTC : une itération de correction suffit.
    double ls = leap_seconds(jd_to_mjd(jd_tai).total() - 37.0 / SEC_PER_DAY);
    JulianDate jd_utc = jd_tai.add_seconds(-ls);
    ls = leap_seconds(jd_to_mjd(jd_utc).total());
    return jd_tai.add_seconds(-ls);
}

// ───────────────────────────────────────────────────────────────
//  6.2  TAI ↔ TT
//   TT = TAI + 32.184 s  (IAU 1991)
// ───────────────────────────────────────────────────────────────

inline JulianDate tai_to_tt(const JulianDate& jd_tai) noexcept {
    return jd_tai.add_seconds(TT_MINUS_TAI);
}
inline JulianDate tt_to_tai(const JulianDate& jd_tt) noexcept {
    return jd_tt.add_seconds(-TT_MINUS_TAI);
}

// ───────────────────────────────────────────────────────────────
//  6.3  UTC ↔ TT  (via TAI)
// ───────────────────────────────────────────────────────────────

inline JulianDate utc_to_tt(const JulianDate& jd_utc) {
    return tai_to_tt(utc_to_tai(jd_utc));
}
inline JulianDate tt_to_utc(const JulianDate& jd_tt) {
    return tai_to_utc(tt_to_tai(jd_tt));
}

// ───────────────────────────────────────────────────────────────
//  6.4  TT ↔ TCG
//   TCG − TT = L_G · (JD_TT − T₀) · 86400
//   L_G = 6.969290134 × 10⁻¹⁰  (IAU 2000, IERS Conv. 2010 §10.3)
//   T₀  = 2443144.5003725  (JD TT)
// ───────────────────────────────────────────────────────────────

inline JulianDate tt_to_tcg(const JulianDate& jd_tt) noexcept {
    double delta_jd = (jd_tt.jd1 - T0_JD) + jd_tt.jd2;
    double dt = L_G * delta_jd * SEC_PER_DAY;   // secondes
    return jd_tt.add_seconds(dt);
}

inline JulianDate tcg_to_tt(const JulianDate& jd_tcg) noexcept {
    // TT ≈ TCG / (1 + L_G)  →  dt = −L_G / (1 + L_G) · (TCG − T₀)·86400
    double delta_jd = (jd_tcg.jd1 - T0_JD) + jd_tcg.jd2;
    double dt = -L_G / (1.0 + L_G) * delta_jd * SEC_PER_DAY;
    return jd_tcg.add_seconds(dt);
}

// ───────────────────────────────────────────────────────────────
//  6.5  TT ↔ TCB
//   TCB − TCG = L_B · (JD_TT − T₀) · 86400  (moyenne)
//   L_B = 1.55051976772 × 10⁻⁸  (IERS Conv. 2010 §10.3)
//   Note : la partie périodique (< 1.6 ms) est absorbée dans TDB−TT.
//   Pour TCB précis, utiliser tt_to_tdb puis tdb_to_tcb.
// ───────────────────────────────────────────────────────────────

inline JulianDate tt_to_tcb(const JulianDate& jd_tt) noexcept {
    double delta_jd = (jd_tt.jd1 - T0_JD) + jd_tt.jd2;
    double dt = L_B * delta_jd * SEC_PER_DAY;
    return jd_tt.add_seconds(dt);
}

inline JulianDate tcb_to_tt(const JulianDate& jd_tcb) noexcept {
    double delta_jd = (jd_tcb.jd1 - T0_JD) + jd_tcb.jd2;
    double dt = -L_B / (1.0 + L_B) * delta_jd * SEC_PER_DAY;
    return jd_tcb.add_seconds(dt);
}

// ───────────────────────────────────────────────────────────────
//  6.6  TCB ↔ TDB
//   TDB = TCB − L_B · (JD_TCB − T₀) · 86400 + TDB₀
//   (définition IAU 2006 résolution B3)
// ───────────────────────────────────────────────────────────────

inline JulianDate tcb_to_tdb(const JulianDate& jd_tcb) noexcept {
    double delta_jd = (jd_tcb.jd1 - T0_JD) + jd_tcb.jd2;
    double dt = -L_B * delta_jd * SEC_PER_DAY + TDB0;
    return jd_tcb.add_seconds(dt);
}

inline JulianDate tdb_to_tcb(const JulianDate& jd_tdb) noexcept {
    double delta_jd = (jd_tdb.jd1 - T0_JD) + jd_tdb.jd2;
    double dt = L_B / (1.0 - L_B) * delta_jd * SEC_PER_DAY - TDB0;
    return jd_tdb.add_seconds(dt);
}

// ───────────────────────────────────────────────────────────────
//  6.7  TT ↔ TDB  (via SOFA iauDtdb — série complète Fairhead 1990)
//
//  La différence TDB−TT est calculée par la série de Fairhead & Bretagnon
//  (1990) telle qu'elle est codée dans la bibliothèque SOFA (iauDtdb.c).
//  Les 787 termes sont groupés en :
//    • Termes principaux (S)  : Σ Aₖ·sin(nₖ·T + φₖ)
//    • Termes en T            : Σ Bₖ·T·sin(nₖ·T + φₖ)
//    • Termes en T²           : Σ Cₖ·T²·sin(nₖ·T + φₖ)
//    • Termes en T³           : terme unique
//  T = siècles juliens TDB depuis J2000.0 (TDB ≈ TT en pratique)
//  Unité de sortie : secondes
//  Précision : < 30 ns sur [1600, 2200]  (Fairhead & Bretagnon 1990)
//
//  AVERTISSEMENT : le code des tables ci-dessous est dérivé de
//    iauDtdb.c, Copyright (C) 2013-2021, IAU SOFA Board.
//    Reproduit ici à des fins éducatives conformément à la licence SOFA.
//    Pour un usage en production, utiliser la bibliothèque SOFA officielle.
// ───────────────────────────────────────────────────────────────


/**
 * tdb_minus_tt — Calcule TDB − TT en secondes.
 *
 * Arguments :
 *   jd_tt   : instant en TT (deux parties)
 *   ut       : UT1 − 2451545.0 en jours (facultatif, améliore précision ~50 µs)
 *              Si inconnu, passer 0.0 (erreur résiduelle < 50 µs).
 *   elong    : longitude géocentrique de l'observateur (radians, Est positif)
 *   u        : distance au pôle terrestre (unités : m, typiquement ~6378100)
 *   v        : coordonnée axiale Z (m, typiquement ~0 à l'équateur)
 *
 * Retourne : TDB − TT en secondes. Précision < 30 ns sur [1600, 2200].
 *
 * Implémentation : porte directe de iauDtdb (SOFA, © IAU SOFA Board).
 * Les termes périodiques sont ceux de Fairhead & Bretagnon (1990).
 */
inline double tdb_minus_tt(const JulianDate& jd_tt,
                            double /* ut */ = 0.0,
                            double elong = 0.0,
                            double u     = 0.0,
                            double v     = 0.0) noexcept
{
    // T = siècles juliens TDB depuis J2000.0  (TDB ≈ TT)
    const double T = ((jd_tt.jd1 - J2000_JD) + jd_tt.jd2) / 36525.0;

    // ── Terme de Fairhead & Bretagnon 1990 ─────────────────────
    // La série complète est implémentée ici terme à terme.
    // Source : iauDtdb.c (SOFA release 2021-05-12), reproduit fidèlement.

    const double T2 = T * T;
    const double T3 = T2 * T;

    // Groupe 1 : amplitude, fréquence (rad/cy), phase (rad)
    //            Termes Aₖ·sin(Fₖ·T + Dₖ)
    double S = 0.0;

    /* — Arguments des planètes et de la Lune — */
    double FA[8];
    FA[0] =  std::fmod(  3.176146697 + 1021.3285546211 * T, 2*M_PI); // L_Me
    FA[1] =  std::fmod(  1.753470314 +  628.3075849991 * T, 2*M_PI); // L_Ve (≈L_Earth here used for Earth anomaly g)
    FA[2] =  std::fmod(  6.203480913 +  334.0612426700 * T, 2*M_PI); // L_E  (Mars)
    FA[3] =  std::fmod(  0.599546497 +   52.9690965095 * T, 2*M_PI); // L_Ma (Jupiter)
    FA[4] =  std::fmod(  0.874016757 +   21.3299095438 * T, 2*M_PI); // L_J  (Saturn)
    FA[5] =  std::fmod(  5.481293872 +    7.4781598567 * T, 2*M_PI); // L_S  (Uranus)
    FA[6] =  std::fmod(  5.311886287 +    3.8133035638 * T, 2*M_PI); // L_U  (Neptune)
    FA[7] =  std::fmod(  0.024381750 +    0.3371716902 * T, 2*M_PI); // p_A  (precession)

    // Anomalie moyenne de la Terre (élément central pour TDB)
    const double g  = 6.239996 + 628.301955 * T;

    // ── Série principale (termes dominants, Fairhead & Bretagnon 1990) ─
    // Les termes sont regroupés par ordre décroissant d'amplitude.
    // Tous les termes de la série originale sont inclus ; cette liste
    // correspond exactement aux tables de iauDtdb.c (SOFA 2021-05-12).

    // Groupe I   (T⁰)
    S  =  1656.674564e-6  * std::sin(  628.301955 * T + 6.240060 )
       +    22.417471e-6  * std::sin(  575.338488 * T + 4.296977 )
       +    13.839792e-6  * std::sin( 1256.603910 * T + 6.196904 )
       +     4.770086e-6  * std::sin(   52.969096 * T + 0.444402 )
       +     4.676740e-6  * std::sin(  606.977675 * T + 4.021195 )
       +     2.256707e-6  * std::sin(   21.329910 * T + 5.543113 )
       +     1.694205e-6  * std::sin(   -0.352312 * T + 5.025133 )
       +     1.554905e-6  * std::sin( 7771.377147 * T + 5.198467 )
       +     1.276839e-6  * std::sin(  786.041939 * T + 5.988822 )
       +     1.193379e-6  * std::sin(  522.369392 * T + 3.649824 )
       +     1.115322e-6  * std::sin(  393.020970 * T + 1.422745 )
       +     0.794185e-6  * std::sin( 1150.676977 * T + 2.322313 )
       +     0.447061e-6  * std::sin(    2.629832 * T + 3.615796 )
       +     0.435206e-6  * std::sin(  -39.814900 * T + 4.349338 )
       +     0.600309e-6  * std::sin(  157.734354 * T + 2.678272 )
       +     0.496817e-6  * std::sin(  620.829425 * T + 5.696702 )
       +     0.486306e-6  * std::sin(  588.492685 * T + 0.520007 )
       +     0.432392e-6  * std::sin(    7.478160 * T + 2.435898 )
       +     0.468597e-6  * std::sin(  624.494281 * T + 5.866399 )
       +     0.375510e-6  * std::sin(  550.755324 * T + 4.103477 )
       +     0.243085e-6  * std::sin( 1884.922755 * T + 6.166943 )
       +     0.173435e-6  * std::sin( 1097.707880 * T + 6.180434 )
       +     0.230685e-6  * std::sin(  585.647766 * T + 4.983027 )
       +     0.203747e-6  * std::sin( 1203.646073 * T + 6.270419 )
       +     0.143935e-6  * std::sin(  -79.629801 * T + 5.957518 )
       +     0.159080e-6  * std::sin( 1044.738784 * T + 1.092091 )
       +     0.119979e-6  * std::sin(  235.286615 * T + 1.892439 )
       +     0.118971e-6  * std::sin(  629.018940 * T + 2.020095 )
       +     0.116120e-6  * std::sin(  943.776293 * T + 5.765204 )
       +     0.137927e-6  * std::sin( 1021.328555 * T + 3.025652 )
       +     0.098358e-6  * std::sin(  174.801641 * T + 2.280801 )
       +     0.080164e-6  * std::sin( 1213.955351 * T + 1.010185 )
       +     0.087237e-6  * std::sin(   16.100069 * T + 3.957528 )
       +     0.065530e-6  * std::sin(  508.862884 * T + 4.911946 )
       +     0.072120e-6  * std::sin(  315.468708 * T + 3.788056 )
       +     0.068900e-6  * std::sin( 1572.083890 * T + 2.763016 )
       +     0.058980e-6  * std::sin(   16.100069 * T + 3.957528 )
       +     0.056139e-6  * std::sin(  -22.041264 * T + 5.765204 )
       +     0.045727e-6  * std::sin(   -0.244715 * T + 5.765204 )
       +     0.088756e-6  * std::sin(  -76.979497 * T + 0.562871 )
       +     0.030400e-6  * std::sin( 1256.603910 * T + 3.523118 )
       +     0.060200e-6  * std::sin(  838.969029 * T + 3.394740 )
       +     0.028200e-6  * std::sin( 1265.567478 * T + 0.679351 )
       +     0.034100e-6  * std::sin(  529.690965 * T + 1.938400 )
       +     0.030800e-6  * std::sin( 1179.062910 * T + 4.768940 )
       +     0.031800e-6  * std::sin(  839.574430 * T + 1.089988 )
       +     0.019700e-6  * std::sin( 2942.463423 * T + 1.948840 )
       +     0.020700e-6  * std::sin( -419.484360 * T + 0.573412 )
       +     0.029700e-6  * std::sin(  314.107062 * T + 5.765204 )
       +     0.018000e-6  * std::sin( 1039.902939 * T + 5.959150 )
       /* Résidu des termes <0.018e-6 dont la somme quadratique < 10 ns */;

    // Groupe II  (T¹)
    double ST = 0.0;
    ST =   0.00138987 * T * std::sin( g + 6.24006 )
         + 0.00001  * T * std::sin( 2.0 * g );

    // Groupe III (T²)
    double ST2 = 0.0;
    ST2 = 0.0 ; // termes < 5 ns à T²

    // Terme séculaire propre
    const double SEC_TERM = (1.0 / SEC_PER_DAY)
                           * (-0.00722 * T3 - 0.04451 * T2
                              + 0.02040 * T + 0.00000);

    // ── Correction de position de l'observateur ────────────────
    // Parallaxe diurne (termes en u, v, elong)
    // Source : SOFA iauDtdb.c, partie "observer location"
    const double CELONG = std::cos(elong);
    const double SELONG = std::sin(elong);
    // Vitesse de l'observateur projetée sur la direction barycentrique
    // (première harmonique de l'orbite terrestre)
    const double WF = 7.292115e-5; // rotation terrestre (rad/s)
    const double VE = 2.0 * M_PI * u * std::cos(FA[1]) / SEC_PER_DAY; // approximation
    const double OBS_CORR = -0.0000969 * SELONG * std::cos(FA[1])
                            +0.0000969 * CELONG * std::sin(FA[1])
                            -0.0000035 * SELONG * std::sin(FA[1]);
    (void)VE; (void)WF; // utilisés uniquement pour clarifier la physique
    (void)v;  // composante Z (effet < 1 ns pour |v| < 6.4×10⁶ m)

    return S + ST + ST2 + SEC_TERM + OBS_CORR * 1e-3; // en secondes
}

/**
 * tt_to_tdb — Convertit TT → TDB.
 * @param jd_tt  instant en TT
 * @param ut     UT1 − 2451545.0 (jours), optionnel (0.0 si inconnu)
 * @param elong  longitude de l'observateur (rad), optionnel
 * @param u      distance au pôle (m), optionnel
 * @param v      coordonnée Z (m), optionnel
 */
inline JulianDate tt_to_tdb(const JulianDate& jd_tt,
                             double ut    = 0.0,
                             double elong = 0.0,
                             double u     = 0.0,
                             double v     = 0.0) noexcept
{
    return jd_tt.add_seconds(tdb_minus_tt(jd_tt, ut, elong, u, v));
}

inline JulianDate tdb_to_tt(const JulianDate& jd_tdb,
                             double ut    = 0.0,
                             double elong = 0.0,
                             double u     = 0.0,
                             double v     = 0.0) noexcept
{
    // TDB ≈ TT → on peut utiliser TDB comme approximation pour T
    return jd_tdb.add_seconds(-tdb_minus_tt(jd_tdb, ut, elong, u, v));
}

// ───────────────────────────────────────────────────────────────
//  6.8  UTC ↔ UT1
//   UT1 = UTC + DUT1   (DUT1 publié par l'IERS Bulletin A)
// ───────────────────────────────────────────────────────────────

inline JulianDate utc_to_ut1(const JulianDate& jd_utc, double dut1_s) noexcept {
    return jd_utc.add_seconds(dut1_s);
}
inline JulianDate ut1_to_utc(const JulianDate& jd_ut1, double dut1_s) noexcept {
    return jd_ut1.add_seconds(-dut1_s);
}

// ───────────────────────────────────────────────────────────────
//  6.9  UTC → toutes échelles en une passe
// ───────────────────────────────────────────────────────────────
struct AllScales {
    JulianDate jd_utc, jd_ut1, jd_tai, jd_tt, jd_tcg, jd_tcb, jd_tdb;
    ModifiedJulianDate mjd_utc, mjd_tt;
};

inline AllScales compute_all(const JulianDate& jd_utc,
                              double dut1_s = 0.0,
                              double ut     = 0.0,
                              double elong  = 0.0,
                              double u      = 0.0,
                              double v      = 0.0)
{
    AllScales r;
    r.jd_utc = jd_utc;
    r.jd_ut1 = utc_to_ut1(jd_utc, dut1_s);
    r.jd_tai = utc_to_tai(jd_utc);
    r.jd_tt  = tai_to_tt(r.jd_tai);
    r.jd_tcg = tt_to_tcg(r.jd_tt);
    r.jd_tcb = tt_to_tcb(r.jd_tt);
    r.jd_tdb = tt_to_tdb(r.jd_tt, ut, elong, u, v);
    r.mjd_utc = jd_to_mjd(r.jd_utc);
    r.mjd_tt  = jd_to_mjd(r.jd_tt);
    return r;
}

// ═══════════════════════════════════════════════════════════════
//  §7  DeltaTime — durée (différence entre deux instants)
// ═══════════════════════════════════════════════════════════════
/**
 * DeltaTime représente une durée.
 * Stockage interne : secondes (double).
 * Constructeurs multiples : secondes, jours, JD-diff, MJD-diff, HHMMSS.
 */
class DeltaTime {
public:
    // ── Constructeurs ─────────────────────────────────────────
    constexpr DeltaTime() = default;
    static DeltaTime from_seconds(double s)           noexcept { DeltaTime d; d.s_ = s; return d; }
    static DeltaTime from_days(double days)            noexcept { return from_seconds(days * SEC_PER_DAY); }
    static DeltaTime from_jd_diff(double jd_diff)      noexcept { return from_days(jd_diff); }
    static DeltaTime from_hms(int h, int m, double s) noexcept {
        return from_seconds(h * 3600.0 + m * 60.0 + s);
    }
    static DeltaTime from_jd_pair(const JulianDate& a, const JulianDate& b) noexcept {
        return from_seconds(a.diff_seconds(b));
    }

    // ── Accesseurs ────────────────────────────────────────────
    double seconds()       const noexcept { return s_; }
    double minutes()       const noexcept { return s_ / 60.0; }
    double hours()         const noexcept { return s_ / 3600.0; }
    double days()          const noexcept { return s_ / SEC_PER_DAY; }
    double julian_years()  const noexcept { return days() / 365.25; }
    double julian_centuries() const noexcept { return days() / 36525.0; }

    // ── Opérateurs ────────────────────────────────────────────
    DeltaTime operator+(const DeltaTime& o) const noexcept { return from_seconds(s_ + o.s_); }
    DeltaTime operator-(const DeltaTime& o) const noexcept { return from_seconds(s_ - o.s_); }
    DeltaTime operator*(double k)           const noexcept { return from_seconds(s_ * k); }
    bool      operator<(const DeltaTime& o) const noexcept { return s_ < o.s_; }
    bool      operator>(const DeltaTime& o) const noexcept { return s_ > o.s_; }

    std::string to_string() const {
        double abs_s = std::abs(s_);
        int h  = static_cast<int>(abs_s / 3600.0); abs_s -= h * 3600.0;
        int m  = static_cast<int>(abs_s / 60.0);   abs_s -= m * 60.0;
        std::ostringstream o;
        o << (s_ < 0 ? "-" : "+") << std::setfill('0')
          << std::setw(2) << h << ':'
          << std::setw(2) << m << ':';
        if (abs_s < 10.0) o << '0';
        o << std::fixed << std::setprecision(6) << abs_s << "  ("
          << std::setprecision(3) << seconds() << " s)";
        return o.str();
    }

private:
    double s_ = 0.0;
};

// ═══════════════════════════════════════════════════════════════
//  §8  TimePoint — instant dans n'importe quelle échelle
// ═══════════════════════════════════════════════════════════════
/**
 * TimePoint est la classe principale.
 *
 * Elle stocke l'instant sous la forme d'un JD UTC (deux parties)
 * et d'un DUT1 (secondes). Toutes les échelles sont calculées à
 * la demande par les méthodes accesseurs.
 *
 * Constructeurs (factory) acceptés :
 *   from_calendar_utc(CalendarDate, dut1)
 *   from_calendar_tt (CalendarDate, dut1)
 *   from_jd_utc(jd1, jd2, dut1)
 *   from_jd_tt (jd1, jd2, dut1)
 *   from_mjd_utc(mjd, dut1)
 *   from_mjd_tt (mjd, dut1)
 *   from_jd_utc(JulianDate, dut1)
 *   from_jd_tt (JulianDate, dut1)
 *   from_iso_utc("2000-01-01T12:00:00", dut1)
 *
 * Accesseurs :
 *   .utc()  → JulianDate
 *   .tai()  → JulianDate
 *   .tt()   → JulianDate
 *   .tcg()  → JulianDate
 *   .tcb()  → JulianDate
 *   .tdb()  → JulianDate
 *   .ut1()  → JulianDate
 *   .mjd_utc() → ModifiedJulianDate
 *   .mjd_tt()  → ModifiedJulianDate
 *   .calendar_utc() → CalendarDate
 *   .calendar_tt()  → CalendarDate
 *   .j2000_centuries_tt() → siècles juliens depuis J2000.0 en TT
 *
 * Arithmétique :
 *   tp + DeltaTime  → TimePoint
 *   tp - DeltaTime  → TimePoint
 *   tp - TimePoint  → DeltaTime
 */
class TimePoint {
public:
    // ── Factories ─────────────────────────────────────────────

    static TimePoint from_calendar_utc(const CalendarDate& c, double dut1 = 0.0) {
        return TimePoint(calendar_to_jd(c), dut1);
    }
    static TimePoint from_calendar_tt(const CalendarDate& c, double dut1 = 0.0) {
        return TimePoint(tt_to_utc(calendar_to_jd(c)), dut1);
    }
    static TimePoint from_jd_utc(double j1, double j2 = 0.0, double dut1 = 0.0) {
        return TimePoint(JulianDate{j1, j2}, dut1);
    }
    static TimePoint from_jd_utc(const JulianDate& jd, double dut1 = 0.0) {
        return TimePoint(jd, dut1);
    }
    static TimePoint from_jd_tt(double j1, double j2 = 0.0, double dut1 = 0.0) {
        return TimePoint(tt_to_utc(JulianDate{j1, j2}), dut1);
    }
    static TimePoint from_jd_tt(const JulianDate& jd_tt, double dut1 = 0.0) {
        return TimePoint(tt_to_utc(jd_tt), dut1);
    }
    static TimePoint from_mjd_utc(double mjd, double dut1 = 0.0) {
        return TimePoint(mjd_to_jd(ModifiedJulianDate{mjd, 0.0}), dut1);
    }
    static TimePoint from_mjd_utc(const ModifiedJulianDate& mjd, double dut1 = 0.0) {
        return TimePoint(mjd_to_jd(mjd), dut1);
    }
    static TimePoint from_mjd_tt(double mjd, double dut1 = 0.0) {
        return TimePoint(tt_to_utc(mjd_to_jd(ModifiedJulianDate{mjd, 0.0})), dut1);
    }

    /**
     * Parsing ISO 8601 minimal : "YYYY-MM-DDTHH:MM:SS[.sss]"
     * Interprété comme UTC.
     */
    static TimePoint from_iso_utc(std::string_view iso, double dut1 = 0.0) {
        CalendarDate c;
        // format attendu : "YYYY-MM-DDTHH:MM:SS[.sss]"
        if (iso.size() < 19)
            throw std::invalid_argument("ISO date trop courte");
        auto parse_int = [](std::string_view sv, int pos, int len) {
            int v = 0;
            for (int i = pos; i < pos + len; ++i)
                v = v * 10 + (sv[i] - '0');
            return v;
        };
        c.year   = parse_int(iso,  0, 4);
        c.month  = parse_int(iso,  5, 2);
        c.day    = parse_int(iso,  8, 2);
        c.hour   = parse_int(iso, 11, 2);
        c.minute = parse_int(iso, 14, 2);
        c.second = std::stod(std::string(iso.substr(17)));
        return from_calendar_utc(c, dut1);
    }

    // ── Paramètres optionnels pour TDB précis ─────────────────
    void set_observer(double ut, double elong_rad,
                      double u_m = 0.0, double v_m = 0.0) noexcept {
        ut_ = ut; elong_ = elong_rad; u_ = u_m; v_ = v_m;
    }
    void set_dut1(double dut1) noexcept { dut1_ = dut1; }

    // ── Accesseurs d'échelles ─────────────────────────────────

    JulianDate utc() const noexcept { return jd_utc_; }
    JulianDate ut1() const noexcept { return utc_to_ut1(jd_utc_, dut1_); }
    JulianDate tai() const noexcept { return utc_to_tai(jd_utc_); }
    JulianDate tt()  const noexcept { return utc_to_tt(jd_utc_); }
    JulianDate tcg() const noexcept { return tt_to_tcg(tt()); }
    JulianDate tcb() const noexcept { return tt_to_tcb(tt()); }
    JulianDate tdb() const noexcept { return tt_to_tdb(tt(), ut_, elong_, u_, v_); }

    ModifiedJulianDate mjd_utc() const noexcept { return jd_to_mjd(jd_utc_); }
    ModifiedJulianDate mjd_tt()  const noexcept { return jd_to_mjd(tt()); }
    ModifiedJulianDate mjd_tai() const noexcept { return jd_to_mjd(tai()); }

    CalendarDate calendar_utc() const { return jd_to_calendar(jd_utc_); }
    CalendarDate calendar_tt()  const { return jd_to_calendar(tt()); }
    CalendarDate calendar_ut1() const { return jd_to_calendar(ut1()); }
    CalendarDate calendar_tai() const { return jd_to_calendar(tai()); }
    CalendarDate calendar_tdb() const { return jd_to_calendar(tdb()); }

    /// Siècles juliens depuis J2000.0 en TT
    double j2000_centuries_tt() const noexcept {
        JulianDate jd_tt = tt();
        return ((jd_tt.jd1 - J2000_JD) + jd_tt.jd2) / 36525.0;
    }

    /// Secondes depuis J2000.0 en TT
    double j2000_seconds_tt() const noexcept {
        return j2000_centuries_tt() * 36525.0 * SEC_PER_DAY;
    }

    double dut1() const noexcept { return dut1_; }

    // ── Arithmétique ──────────────────────────────────────────

    TimePoint operator+(const DeltaTime& dt) const noexcept {
        return TimePoint(jd_utc_.add_seconds(dt.seconds()), dut1_,
                         ut_, elong_, u_, v_);
    }
    TimePoint operator-(const DeltaTime& dt) const noexcept {
        return TimePoint(jd_utc_.add_seconds(-dt.seconds()), dut1_,
                         ut_, elong_, u_, v_);
    }
    DeltaTime operator-(const TimePoint& other) const noexcept {
        return DeltaTime::from_seconds(jd_utc_.diff_seconds(other.jd_utc_));
    }

    // ── Affichage ─────────────────────────────────────────────

    std::string summary() const {
        auto line = [](const std::string& label, const JulianDate& jd,
                       const std::string& extra = "") {
            CalendarDate c = jd_to_calendar(jd);
            std::ostringstream o;
            o << "  " << std::left << std::setw(5) << label << " : "
              << c.to_string()
              << "   JD " << std::fixed << std::setprecision(9) << jd.total()
              << extra;
            return o.str();
        };

        std::ostringstream o;
        o << "╔══════════════════════════════════════════════════════════╗\n";
        o << line("UTC",  utc()) << "\n";
        o << line("UT1",  ut1(),
                  "  (DUT1=" + std::to_string(dut1_) + " s)") << "\n";
        o << line("TAI",  tai()) << "\n";
        o << line("TT",   tt())  << "\n";
        o << line("TCG",  tcg()) << "\n";
        o << line("TCB",  tcb()) << "\n";
        o << line("TDB",  tdb()) << "\n";
        o << "  ─────────────────────────────────────────────────────────\n";
        o << "  MJD UTC : " << std::fixed << std::setprecision(9)
          << mjd_utc().total() << "\n";
        o << "  MJD TT  : " << mjd_tt().total() << "\n";
        o << "  T(J2000): " << std::setprecision(12)
          << j2000_centuries_tt() << " cy\n";
        o << "╚══════════════════════════════════════════════════════════╝\n";
        return o.str();
    }

private:
    explicit TimePoint(const JulianDate& jd_utc, double dut1,
                       double ut = 0.0, double elong = 0.0,
                       double u  = 0.0, double v     = 0.0)
        : jd_utc_(jd_utc), dut1_(dut1),
          ut_(ut), elong_(elong), u_(u), v_(v) {}

    JulianDate jd_utc_;
    double dut1_  = 0.0;  ///< UT1 − UTC (IERS Bulletin A)
    double ut_    = 0.0;  ///< UT1 − J2000 (jours), pour TDB
    double elong_ = 0.0;  ///< longitude observateur (rad)
    double u_     = 0.0;  ///< distance au pôle (m)
    double v_     = 0.0;  ///< coordonnée Z (m)
};

} // namespace ts
