MRF.correctBoundaryVelocity(U);

// Calculate the drag force
volScalarField KInv
(
    (scalar(180)/sqr(DAS))*sqr(min(alpha1, 0.95))/pow3(max(alpha2, 0.05))
);

volScalarField Kd("Kd", rho*mu2/rho2*KInv);

// Update the kinematic density
rhok = rho2*(scalar(1) - betaT*(T - TRef) - betaC*(C2 - CRef));


fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(rho*mu2/rho2, U)
 ==
    fvOptions(rho, U)
  - fvm::Sp(Kd, U)
);
UEqn.relax();
fvOptions.constrain(UEqn);


if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
