{
    // Artificial limits for eutectic reaction
    volScalarField condition = pos0(C - Calpha);

    volScalarField Thigh = Tsol + condition*(Teut + Tiso - Tsol);

    volScalarField Tlow  = Tsol + condition*(Teut - Tiso - Tsol);

    // All liquid
    condition = pos0(T - Tliq);
    {
        alpha10 *= scalar(1) - condition;

        dFdT *= scalar(1) - condition;

        T0 = T0 + condition*(T - T0);
    }

    // Primary solidification
    condition = pos(Tliq - T) * pos(T - Thigh);
    {
        visitedEut *= scalar(1) - condition;

        dFdT = dFdT
        + condition *
        (
            (Tliq - Tmelt) / (sqr(T - Tmelt)) / (scalar(1) - kp)
        - dFdT
        );

        T0 = T0
        + condition *
        (
            (Tliq - alpha1*(scalar(1) - kp)*Tmelt)
          / (scalar(1) - alpha1*(scalar(1) - kp))
        - T0
        );
    }

    // Possible eutectic solidification (linear regression)
    {
        condition.primitiveFieldRef() =
            pos0(Thigh - T)*pos(T - Tlow)*neg(visitedEut - 0.5);
        {
            alphaEut = alphaEut + condition*(alpha1 - alphaEut);

            visitedEut = visitedEut + condition*(scalar(1) - visitedEut);
        }

        condition = pos0(Thigh - T ) * pos(T - Tlow);

        dFdT.primitiveFieldRef() = dFdT
        + condition *
        (
            (alphaEut - scalar(1))
           / stabilise(Thigh - Tlow, dimensionedScalar(dimTemperature, small))
        - dFdT
        );

        T0 = T0
        + condition *
        (
            Tlow
          + (alpha1 - scalar(1))
          / stabilise(dFdT, dimensionedScalar(dimless/dimTemperature, small))
        - T0
        );
    }

    // All solid
    condition = pos0(Tlow - T);
    {
        alpha10 = alpha10 + condition*(scalar(1) - alpha10);

        dFdT *= scalar(1) - condition;

        T0 = T0 + condition*(T - T0);
    }
}

dFdT.correctBoundaryConditions();
T0.boundaryFieldRef() = T.boundaryField();
