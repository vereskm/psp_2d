function power = validateModePower(complexPower,context)
%VALIDATEMODEPOWER Validate modal power before normalizing a port mode.

arguments
    complexPower (1,1) double
    context (1,1) string = "port mode"
end

assert(isfinite(real(complexPower)) && isfinite(imag(complexPower)), ...
    "%s power is not finite.",context);

power = real(complexPower);
assert(power > 0, ...
    "%s has non-positive forward power: %.16g.",context,power);

relativeImaginaryPower = abs(imag(complexPower))/power;
assert(relativeImaginaryPower <= 1e-8, ...
    "%s power is significantly complex (imaginary/real = %.3g). " + ...
    "Check the port mode, PML overlap, and propagation direction.", ...
    context,relativeImaginaryPower);
end
