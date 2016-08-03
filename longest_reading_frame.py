#!/bin/env python
import sys

def longest_frame():
    max_frame_length = -1
    current_frame_length = 0
    in_open_frame = False
    position = 0
    max_position = 0

    for line in sys.stdin:
        line = line.strip()
        for char in line:
            if char == "M":
                current_frame_length = 1
                in_open_frame = True
            elif char == "-":
                if in_open_frame and current_frame_length > max_frame_length:
                    max_frame_length = current_frame_length
                    max_position = position - max_frame_length
                current_frame_length = 0
                in_open_frame = False
            elif in_open_frame:
                current_frame_length += 1

            position += 1

    print "%s (%s of %s)" % (max_frame_length, max_position, position)


if __name__ == "__main__":
    longest_frame()
