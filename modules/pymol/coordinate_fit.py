"""Least-squares rigid fit of paired Cartesian coordinates."""

import numpy as np


def rigid_fit(mobile, reference):
    mobile = np.asarray(mobile, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    moving = mobile - mobile_center
    target = reference - reference_center
    left, _, right = np.linalg.svd(moving.T @ target)
    if np.linalg.det(left @ right) < 0:
        left[:, -1] *= -1
    rotation = left @ right
    translation = reference_center - mobile_center @ rotation
    residual = moving @ rotation - target
    rmsd = float(np.sqrt(np.mean(np.sum(residual * residual, axis=1))))
    return rotation, translation, rmsd
