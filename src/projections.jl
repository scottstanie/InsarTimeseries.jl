
function convert_xyz_latlon_to_enu(lat_lons::Array{<:AbstractFloat, 2}, xyz_array::Array{<:AbstractFloat, 2})
    out = similar(xyz_array)
    for i = 1:size(lat_lons, 2)
        out[:, i] = rotate_xyz_to_enu(xyz_array[:, i], lat_lons[:, i]...)
    end
    return out
end


"""Rotates a vector in XYZ coords to ENU

Args:
    xyz: 3xk array of xyz vectors
    lat: latitude (deg) of point to rotate into
    lon: longitude (deg) of point to rotate into

Reference: https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates

Example:
>>> rotate_xyz_to_enu([-2, -3, 1], 0, 0)
[ 1.,  2.,  3.]
"""
function rotate_xyz_to_enu(xyz, lat, lon)
    # Rotate about axis 3 with longitude, then axis 1 with latitude
    R3 = rot(90 + lon, 3, in_degrees=true)
    R1 = rot(90 - lat, 1, in_degrees=true)
    return (R3 * R1) * xyz
end



"""
Find a 3x3 euler rotation matrix given an angle and axis.

Rotation matrix used for rotating a vector about a single axis.

Args:
    `angle` in degrees to rotate
    axis: 1, 2 or 3 for R1, R2, or R3
    in_degrees (bool): specify the angle in degrees. if false, using radians
"""
function rot(angle, axis; in_degrees=true)
    R = Matrix{Float64}(I, 3, 3)
    if in_degrees
        angle = deg2rad(angle)
    end
    cang = cos(angle)
    sang = sin(angle)
    if axis == 1
        R[2, 2] = cang
        R[3, 3] = cang
        R[2, 3] = sang
        R[3, 2] = -sang
    elseif axis == 2
        R[1, 1] = cang
        R[3, 3] = cang
        R[1, 3] = -sang
        R[3, 1] = sang
    elseif axis == 3
        R[1, 1] = cang
        R[2, 2] = cang
        R[2, 1] = -sang
        R[1, 2] = sang
    else
        throw("axis must be 1, 2 or 2")
    end
    return R
end

"""Find magnitude of an ENU vector in the LOS direction

Rotates the line of sight vector to ENU coordinates at
(lat, lon), then dots with the enu data vector

Args:
    enu (list[float], ndarray[float]): E,N,U coordinates, either
        as list of 3, or a (3, k) array of k ENU vectors
    los_vec (ndarray[float]) length 3 line of sight, in XYZ frame
    lat (float): degrees latitude of los point
    lon (float): degrees longitude of los point
    enu_coeffs (ndarray) size 3 array of the E,N,U coefficients
    of a line of sight vector. Comes from `find_enu_coeffs`.
        If this arg is used, others are not needed

Returns:
    ndarray: magnitudes same length as enu input, (k, 1)

Examples:
>>> print('%.2f' % project_enu_to_los([1,2,3],[1, 0, 0], 0, 0))
-2.00
>>> print('%.2f' % project_enu_to_los([1,2,3],[0, 1, 0], 0, 0))
-3.00
>>> print('%.2f' % project_enu_to_los([1,2,3],[0, 0, 1], 0, 0))
1.00
"""
function project_enu_to_los(enu, los_vec=nothing, lat=nothing, lon=nothing, enu_coeffs=nothing)
    if isnothing(enu_coeffs)
        los_hat = los_vec / norm(los_vec)
        enu_coeffs = rotate_xyz_to_enu(los_hat, lat, lon)
    end
    return enu_coeffs * enu
end
