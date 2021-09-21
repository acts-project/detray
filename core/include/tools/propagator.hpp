/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <iostream>

namespace detray {

/** A sterily track inspector instance */
struct void_track_inspector {

  /** Void operator */
  template <typename track_t> void operator()(const track_t & /*ignored*/) {
    return;
  }
};

/** Tempalted propagator class, using a
 *
 * @tparam stepper_t for the transport
 * @tparam navigator_t for the navigation
 *
 **/
template <typename stepper_t, typename navigator_t> struct propagator {

  stepper_t _stepper;
  navigator_t _navigator;

  /** Can not be default constructed */
  propagator() = delete;
  /** Only valid constructor with a
   * @param s stepper
   * @param n navigator
   * by move semantics
   **/
  propagator(stepper_t &&s, navigator_t &&n)
      : _stepper(std::move(s)), _navigator(std::move(n)) {}

  /** Propagate method
   *
   * @tparam track_t is the type of the track
   *
   * @param t_in the track at input
   *
   * @return a track at output
   */
  template <typename track_t, typename track_inspector_t>
  track_t propagate(const track_t &t_in, track_inspector_t &t_inspector) const {

    track_t t_out(t_in);
    typename stepper_t::state s_state(t_out);
    typename navigator_t::state n_state;

    bool heartbeat = _navigator.status(n_state, s_state());
    t_inspector(t_out);
    // Run while there is a heartbeat
    while (heartbeat) {
      // (Re-)target
      heartbeat &= _navigator.target(n_state, s_state());
      // Take the step
      heartbeat &= _stepper.step(s_state, n_state());
      t_inspector(t_out);
      // And check the status
      heartbeat &= _navigator.status(n_state, s_state());
    }
    return t_out;
  };
};

} // namespace detray
