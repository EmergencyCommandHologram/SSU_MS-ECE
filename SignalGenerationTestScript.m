% Assume vectors
I1 = out.CH1_I;
Q1 = out.CH1_Q;
I_Dop = out.Dop_I;
Q_Dop = out.Dop_Q;
I_ref = out.CH1_IRef;
Q_ref = out.CH1_QRef;

% Received signal (delayed chirp)
rx = I1 + 1j*Q1;
ref = I_ref - 1j*Q_ref;

beat = rx .* conj(ref);

phase = unwrap(angle(beat));
plot(real(beat))

