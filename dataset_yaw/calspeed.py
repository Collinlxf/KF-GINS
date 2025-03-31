import math

speed = 4.61 # m/s
heading = 319.52 * math.pi / 180.0  # rad

vn = speed * math.cos(heading)
ve = speed * math.sin(heading)

print(f"北向速度 (Vn): {vn} m/s")
print(f"东向速度 (Ve): {ve} m/s")
