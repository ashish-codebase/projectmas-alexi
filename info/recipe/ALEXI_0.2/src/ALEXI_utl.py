def get_pbl_table(profile, z1=50):
    # find the temperature at z1
    heights = np.array(
        [0.0, 100.0, 300.0, 500.0, 700.0, 1000.0, 1400.0, 1800.0, 2200.0, 2600.0, 3000.0, 3500.0, 4000.0, 4500.0])
    temp_z1 = np.interp(z1, heights, profile)
    new_temp_profile = profile[np.where(heights > z1)[0]]
    new_temp_profile = np.insert(new_temp_profile, 0, temp_z1)

    new_heights = heights[np.where(heights > z1)[0]]
    new_heights = np.insert(new_heights, 0, z1)

    dtheta = 0.05
    tlast = new_temp_profile[0]
    i = 0
    zlast = 50.
    tint = 0
    tabtheta = []
    tabz2 = []
    for it in range(8000):
        t = tlast + dtheta
        while t > new_temp_profile[i + 1]:
            i += 1
            if (i > 12):
                break
        if (i > 12):
            break
        f = (t - new_temp_profile[i]) / (new_temp_profile[i + 1] - new_temp_profile[i])
        z = new_heights[i] + f * (new_heights[i + 1] - new_heights[i])
        dz = z - zlast
        tint = tint + 0.5 * dz * (t + tlast)  # Trapezoid rule integration
        tabtheta.append(tint)
        tabz2.append(z)
        tlast = t
        zlast = z

    return tabtheta, tabz2


def get_geopotential_height(press_level):
    return 8 * (1014. - press_level)


def get_theta_sum(profile, th2, z1=50):
    # find the temperature at z1
    heights = np.array(
        [0.0, 100.0, 300.0, 500.0, 700.0, 1000.0, 1400.0, 1800.0, 2200.0, 2600.0, 3000.0, 3500.0, 4000.0, 4500.0])
    temp_z1 = np.interp(z1, heights, profile)
    new_temp_profile = profile[np.where(heights > z1)[0]]
    new_temp_profile = np.insert(new_temp_profile, 0, temp_z1)

    new_heights = heights[np.where(heights > z1)[0]]
    new_heights = np.insert(new_heights, 0, z1)

    prof_cdf = integrate.cumtrapz(new_temp_profile, new_heights)
    prof_cdf = np.insert(prof_cdf, 0, 0)
    z2 = np.interp(th2, profile, heights)
    thetasum = np.interp(z2, new_heights, prof_cdf)
    
    return thetasum, temp_z1, z2
