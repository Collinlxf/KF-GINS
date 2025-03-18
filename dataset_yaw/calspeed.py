import math

speed = 22.59  # m/s
heading = 3.08 * math.pi / 180.0  # rad

vn = speed * math.cos(heading)
ve = speed * math.sin(heading)

print(f"北向速度 (Vn): {vn} m/s")
print(f"东向速度 (Ve): {ve} m/s")
