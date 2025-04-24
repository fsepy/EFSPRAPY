import numpy as np


# Assume solid_angle_triangle(P, A, B, C) function is defined as before

def solid_angle_triangle(P, A, B, C):
    """
    Calculates the solid angle subtended by triangle ABC at point P.
    Uses the Oosterom and Strackee formula (simplified).
    Note: Be cautious with numerical stability and sign conventions.
    """
    vec_a = A - P
    vec_b = B - P
    vec_c = C - P

    mag_a = np.linalg.norm(vec_a)
    mag_b = np.linalg.norm(vec_b)
    mag_c = np.linalg.norm(vec_c)

    # Normalize vectors (optional but simplifies formula slightly if done here)
    # norm_a = vec_a / mag_a
    # norm_b = vec_b / mag_b
    # norm_c = vec_c / mag_c
    # Then use these in the formula below with mags = 1

    # Triple product numerator term
    numerator = np.abs(np.dot(np.cross(vec_a, vec_b), vec_c))

    # Denominator terms (using original vectors and magnitudes)
    denom = (mag_a * mag_b * mag_c +
             np.dot(vec_a, vec_b) * mag_c +
             np.dot(vec_b, vec_c) * mag_a +
             np.dot(vec_c, vec_a) * mag_b)

    # Ensure denominator is not zero or problematic
    if np.isclose(denom, 0):
        # Handle degenerate cases, e.g., P lies on the plane of ABC
        # Or if P is very far, the angle might be close to zero
        # A more robust check is needed here.
        # For simplicity, return 0, but this might not be correct
        # if P is inside the triangle projection or on an edge.
        print("Warning: Denominator close to zero.")
        # A better approach might be to calculate angles directly if denom is small
        # For now, just return zero angle if numerically unstable.
        # This calculation also breaks if P coincides with A, B, or C.
        if np.isclose(mag_a * mag_b * mag_c, 0):  # Check if P coincides with a vertex
            return 0.0  # Or handle appropriately

        # If numerator is also zero, likely P is on the plane OUTSIDE triangle.
        if np.isclose(numerator, 0):
            return 0.0
        else:
            # If P is on the plane INSIDE the triangle, solid angle should be 2*pi
            # This needs a point-in-triangle test in the plane. Not implemented here.
            # For now, return NaN or raise error for unhandled case.
            return np.nan  # Indicates potential issue

    # Calculate the angle using atan2 for better quadrant handling
    # Using atan2(numerator_signed, denom) is often more robust
    # Need the signed triple product for atan2:
    triple_product_signed = np.dot(np.cross(vec_a, vec_b), vec_c)
    omega = 2 * np.arctan2(triple_product_signed, denom)

    # Ensure omega is positive (solid angle is non-negative)
    # The formula might give results between -2pi and 2pi.
    # A positive solid angle is expected.
    # If P is 'behind' the triangle relative to its normal, the sign might flip.
    # Often |omega| is taken, but check reference for precise definition.
    return np.abs(omega)  # Often the magnitude is desired


def view_factor_point_to_polygon(P, vertices):
    """
    Calculates the view factor from point P to a polygon defined by vertices.
    Vertices should be ordered (clockwise or counter-clockwise).
    Uses fan triangulation from the first vertex.
    """
    num_vertices = len(vertices)
    if num_vertices < 3:
        return 0.0  # Not a polygon

    total_solid_angle = 0.0
    V1 = vertices[0]

    for i in range(1, num_vertices - 1):
        V_i = vertices[i]
        V_i_plus_1 = vertices[i + 1]

        # Calculate solid angle for triangle (V1, V_i, V_i+1)
        omega_triangle = solid_angle_triangle(P, V1, V_i, V_i_plus_1)

        if np.isnan(omega_triangle):
            print(f"Warning: NaN encountered for triangle P-{V1}-{V_i}-{V_i_plus_1}. Skipping.")
            # Or handle error more robustly
            continue  # Or return NaN/error for the whole polygon

        total_solid_angle += omega_triangle

    # Ensure total solid angle is physically meaningful (e.g., <= 4*pi)
    # Clamping might hide issues, but can prevent nonsensical results.
    # total_solid_angle = np.clip(total_solid_angle, 0, 4 * np.pi)

    F_P_to_Polygon = total_solid_angle / (4 * np.pi)
    return F_P_to_Polygon


if __name__ == '__main__':
    # --- Example Usage ---
    P = np.array([0.0, 0.0, 0.0])

    # Define a square on the z=1 plane
    vertices_square = [
        np.array([2.0, -1.0,    11.0]),  # V1
        np.array([2.0, 1.0,     11.0]),  # V2
        np.array([-2.0, 1.0,    11.0]),  # V3
        np.array([-2.0, -1.0,   11.0])  # V4
    ]

    F_P_to_Square = view_factor_point_to_polygon(P, vertices_square)
    print(f"View Factor F_(P -> Square): {F_P_to_Square*4:.6f}")

    # Define a more complex polygon (ensure correct ordering)
    # vertices_poly = [
    #     np.array([2.0, 0.0, 1.0]),
    #     np.array([1.0, 2.0, 1.5]),
    #     np.array([-1.0, 1.0, 1.0]),
    #     np.array([-1.0, -1.0, 1.2]),
    #     np.array([0.0, -2.0, 0.8])
    # ]
    #
    # F_P_to_Poly = view_factor_point_to_polygon(P, vertices_poly)
    # print(f"View Factor F_(P -> Polygon): {F_P_to_Poly:.6f}")
