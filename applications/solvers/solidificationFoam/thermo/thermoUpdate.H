{
    // Update the phase fraction due to heat/solute transport
    alpha1 = alpha10 + dFdT*(T - T0);

    alpha1 = min(max(alpha1, scalar(0)), scalar(1));

    // Update the phase compositions
    C2 = C / (scalar(1) + alpha1*(kp-scalar(1)));

    C2 = min(C2, Ceut);

    C1 = kp*C2;
}

alpha1.correctBoundaryConditions();

alpha2 = scalar(1) - alpha1;

C1.correctBoundaryConditions();
C2.correctBoundaryConditions();
