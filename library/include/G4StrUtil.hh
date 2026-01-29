/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4StrUtil.hh
/// \brief Namespace defined only for Geant4 v10, which provides calls to
///	G4String utililty functions.  This allows application code to be
///	migrated to Geant4 v11, and remain compatible with Geant4 v10.
//
//  20251116  G4CMP-525 -- Adapted from SuperCDMS version.
//  20251203  G4CMP-551 -- Modify implementations of *strip() to match G4String

#ifndef G4StrUtil_hh
#define G4StrUtil_hh 1

#include "G4Version.hh"

#if G4VERSION_NUMBER < 1100
#include "G4String.hh"

namespace G4StrUtil {
  inline void to_lower(G4String& str) { str.toLower(); }

  inline G4String to_lower_copy(G4String str) { str.toLower(); return str; }

  inline void to_upper(G4String& str) { str.toUpper(); }

  inline G4String to_upper_copy(G4String str) { str.toUpper(); return str; }

  inline G4String lstrip_copy(G4String str, char ch=' ') {
    return str.strip(G4String::leading, ch);
  }

  inline G4String rstrip_copy(G4String str, char ch=' ') {
    return str.strip(G4String::trailing, ch);
  }

  inline G4String strip_copy(G4String str, char ch=' ') {
    return str.strip(G4String::both, ch);
  }

  inline void lstrip(G4String& str, char ch=' ') {
    str = lstrip_copy(str, ch);
  }

  inline void rstrip(G4String& str, char ch=' ') {
    str = rstrip_copy(str, ch);
  }

  inline void strip(G4String& str, char ch=' ') {
    str = strip_copy(str, ch);
  }

  inline G4bool contains(const G4String& str, const G4String& ss) {
    return str.contains(ss);
  }

  inline G4bool contains(const G4String& str, char ss) {
    return str.contains(ss);
  }

  inline G4bool contains(const G4String& str, const char* ss) {
    return str.contains(ss);
  }

  inline G4int icompare(const G4String& lhs, const G4String& rhs) {
    return lhs.compareTo(rhs, G4String::ignoreCase);
  }

  inline std::istream& readline(std::istream& is, G4String& str,
                                G4bool skipWhite=true) {
    return str.readLine(is, skipWhite);
  }
}

#endif  /* G4 10 */
#endif  /* G4StrUtil_hh */
