/*
* Brian R Taylor
* brian.taylor@bolderflight.com
* 
* Copyright (c) 2022 Bolder Flight Systems Inc
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the “Software”), to
* deal in the Software without restriction, including without limitation the
* rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
* sell copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
* IN THE SOFTWARE.
*/

#ifndef UNITS_SRC_CONVVEL_H_  // NOLINT
#define UNITS_SRC_CONVVEL_H_

#if defined(ARDUINO)
#include <Arduino.h>
#endif
#include <type_traits>

namespace bfs {
/* Units for measuring linear velocity */
enum class LinVelUnit {
  FPS,  // feet per second, ft/s
  MPS,  // meters per second, m/s
  KPS,  // kilometers per second, km/s
  IPS,  // inches per second. in/s
  KPH,  // kilometers per hour, km/h
  MPH,  // miles per hour, mi/h
  KTS,  // knots
  FPM   // feet per minute, ft/min
};
/* 
* Utility to convert between linear velocity units:
* Input the value to convert, the unit the value is currently in, and the unit
* you are converting to, i.e. 'convvel(1, LinVelUnit::FPS, LinVelUnit::MPS)'
* converts 1 ft/s to m/s.
*/
template<typename T>
T convvel(const T val, const LinVelUnit input, const LinVelUnit output) {
  static_assert(std::is_floating_point<T>::value,
              "Only floating point types supported");
  /* Trivial case where input and output units are the same */
  if (input == output) {return val;}
  /* Convert input to SI */
  T in_val;
  switch (input) {
    case LinVelUnit::FPS: {
      in_val = val * static_cast<T>(0.3048);
      break;
    }
    case LinVelUnit::MPS: {
      in_val = val;
      break;
    }
    case LinVelUnit::KPS: {
      in_val = val * static_cast<T>(1000);
      break;
    }
    case LinVelUnit::IPS: {
      in_val = val * static_cast<T>(0.0254);
      break;
    }
    case LinVelUnit::KPH: {
      in_val = val * static_cast<T>(1000) / static_cast<T>(3600);
      break;
    }
    case LinVelUnit::MPH: {
      in_val = val * static_cast<T>(1609.344) / static_cast<T>(3600);
      break;
    }
    case LinVelUnit::KTS: {
      in_val = val * static_cast<T>(1852) / static_cast<T>(3600);
      break;
    }
    case LinVelUnit::FPM: {
      in_val = val * static_cast<T>(0.3048) / static_cast<T>(60);
      break;
    }
  }
  /* Convert to output */
  T out_val;
  switch (output) {
    case LinVelUnit::FPS: {
      out_val = in_val / static_cast<T>(0.3048);
      break;
    }
    case LinVelUnit::MPS: {
      out_val = in_val;
      break;
    }
    case LinVelUnit::KPS: {
      out_val = in_val / static_cast<T>(1000);
      break;
    }
    case LinVelUnit::IPS: {
      out_val = in_val / static_cast<T>(0.0254);
      break;
    }
    case LinVelUnit::KPH: {
      out_val = in_val / static_cast<T>(1000) * static_cast<T>(3600);
      break;
    }
    case LinVelUnit::MPH: {
      out_val = in_val / static_cast<T>(1609.344) * static_cast<T>(3600);
      break;
    }
    case LinVelUnit::KTS: {
      out_val = in_val / static_cast<T>(1852) * static_cast<T>(3600);
      break;
    }
    case LinVelUnit::FPM: {
      out_val = in_val / static_cast<T>(0.3048) * static_cast<T>(60);
      break;
    }
  }
  return out_val;
}

}  // namespace bfs

#endif  // UNITS_SRC_CONVVEL_H_ NOLINT
