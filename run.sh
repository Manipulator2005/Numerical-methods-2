#!/bin/bash

set -e

meson setup build-release --buildtype=release && 
meson compile -C build-release && 
ln -sf build-release/prog prog &&
./prog