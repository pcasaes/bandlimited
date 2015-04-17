Generates a bandlimited waveform (square, pulse, saw, reverse saw, triangle and saw-triangle). It's purpose is to produce a signal without aliasing.

bandlimited~ works by using a series of wavetables with different quanities of harmonics. This is done to keep CPU usage at a minimum. The wavetable with the highest harmonic content has a maximum of 1104 harmonics. Any frequency below 20hz at 44.1kHz will start to geometrically use more CPU power.
