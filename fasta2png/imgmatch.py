import cv2
import numpy as np
import matplotlib.pyplot as plt

# example fourier transform
img = cv2.imread('NC_008253_h_seq2.png', 0)
f = np.fft.fft2(img)
fshift = np.fft.fftshift(f)
magnitude_spectrum = 20*np.log(np.abs(fshift))
cv2.imwrite("mag_spec_host.png")

# img read
img_vir_mag = cv2.imread("mag_spec.png")
img_vir_mag_gray = cv2.cvtColor(img_vir_mag, cv2.COLOR_BGR2GRAY)
img_host_mag = cv2.imread("mag_spec_host.png")
img_host_mag_gray = cv2.cvtColor(img_host_mag, cv2.COLOR_BGR2GRAY)

# detection
orb = cv2.ORB_create()
kp_v, des_v = orb.detectAndCompute(img_vir_mag, None)
kp_h, des_h = orb.detectAndCompute(img_host_mag, None)
# example keypoint drawing
virus_keypoints_draw = cv2.drawKeypoints(img_vir_mag_gray, kp_v, img_vir_mag)

# matching
bf_hamm = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
matches_hamm = bf_hamm.match(des_v, des_h)
matches_hamm = sorted(matches_hamm, key=lambda x:x.distance)
match_hamm = cv2.drawMatches(img_vir_mag, kp_v, img_host_mag, kp_h, matches_hamm[:10], None, flags=cv2.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS)

# show matches
plt.imshow(match_hamm),plt.show()
