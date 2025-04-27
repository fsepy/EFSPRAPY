from fsetools.lib.fse_thermal_radiation import phi_parallel_any_br187


def goal_seek_binary(
        func, target_output, initial_guess=1.0, step_size=0.1, max_steps=1000, tolerance=1e-6, max_binary_iterations=50
):
    """
    Find the input value that produces a target output using incremental steps
    followed by binary search. Only searches for values greater than or equal to initial_guess.
    Designed to work with inversely proportional functions where output decreases as input increases.

    Args:
        func: The function to goal seek (takes one argument, returns one value)
        target_output: The desired output value
        initial_guess: Starting point for the search (minimum allowed value)
        step_size: Initial step size for finding bounds
        max_steps: Maximum number of steps to find bounds
        tolerance: Acceptable difference between actual and target output
        max_binary_iterations: Maximum number of binary search iterations

    Returns:
        The input value that produces the target output (or closest approximation)
    """
    # Start with the initial guess
    x = initial_guess
    initial_output = func(x)

    # Initialize bound variables
    lower_bound = x  # Set lower bound to initial_guess to enforce constraint
    upper_bound = None
    lower_output = initial_output
    upper_output = None

    # Check if initial guess is already close enough
    if abs(initial_output - target_output) <= tolerance:
        return x

    # For inverse proportional functions, we need to check which direction to move
    # If target is less than current output, we need to increase input
    # If target is greater than current output, we can't reach it by increasing input
    if target_output > initial_output:
        raise ValueError(f"Target output {target_output} cannot be reached with values >= {initial_guess} "
                         f"for an inversely proportional function. "
                         f"Initial output at {initial_guess} is {initial_output}.")

    # Phase 1: Find upper bound where output is lower than target
    steps_taken = 0
    x = initial_guess  # Reset x to initial_guess

    while steps_taken < max_steps:
        # Take a step
        x += step_size
        output = func(x)

        # For inverse proportional functions, we're looking for output <= target
        if output <= target_output:
            upper_bound = x
            upper_output = output
            break
        else:
            # Update lower bound and continue
            lower_bound = x
            lower_output = output
            # Accelerate the step size for faster convergence
            step_size *= 2

        steps_taken += 1

    # Check if we found an upper bound
    if upper_bound is None:
        raise ValueError(f"Could not find an input value >= {initial_guess} that produces an output <= {target_output} "
                         f"within {max_steps} steps.")

    # Phase 2: Binary search between the bounds
    binary_iterations = 0

    while binary_iterations < max_binary_iterations:
        # Check if we're close enough
        if abs(upper_bound - lower_bound) < tolerance or abs(upper_output - target_output) < tolerance or abs(
                lower_output - target_output) < tolerance:
            # Return the closer bound
            if abs(upper_output - target_output) < abs(lower_output - target_output):
                return upper_bound
            else:
                return lower_bound

        # Calculate midpoint
        mid = (lower_bound + upper_bound) / 2
        mid_output = func(mid)

        # Check if midpoint is close enough
        if abs(mid_output - target_output) < tolerance:
            return mid

        # Update bounds - note the reversed logic for inverse proportional functions
        if mid_output > target_output:
            lower_bound = mid
            lower_output = mid_output
        else:
            upper_bound = mid
            upper_output = mid_output

        binary_iterations += 1

    # Return the closer approximation after max iterations
    if abs(upper_output - target_output) < abs(lower_output - target_output):
        return upper_bound
    else:
        return lower_bound


def my_function(x):
    return phi_parallel_any_br187(9, 3, 9 / 2, 3 / 2, x)


def sep_parallel_any_br187(w, h, q_emit, q_crit):
    phi_target = q_crit / q_emit
    result = goal_seek_binary(
        lambda sep: phi_parallel_any_br187(w, h, w / 2, h / 2, sep),
        phi_target,
        initial_guess=0.001
    )
    return result
