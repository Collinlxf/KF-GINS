import math

speed = 6.38  # m/s
heading = 337.33 * math.pi / 180.0  # rad

vn = speed * math.cos(heading)
ve = speed * math.sin(heading)

print(f"北向速度 (Vn): {vn} m/s")
print(f"东向速度 (Ve): {ve} m/s")
