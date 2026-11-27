#!/usr/bin/env bash
# One-time setup of everything the decoders need, on a fresh Ubuntu 24.04 PC (x86-64).
#
#   cd ~/fbs-upload && bash setup_ubuntu.sh
#
# Installs: build tools and Python packages, libbladeRF + bladeRF-cli (built from source into /usr/local, as on
# the Pi) with the BladeRF 2.0 FPGA images, GNU Radio + gr-gsm (GSM), and srsRAN 4G with our patch (LTE).
# Takes 15-30 minutes. Safe to run again (it skips what is already done). Asks for your password once (sudo).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="$(nproc)"
step() { printf '\n\033[1;36m== %s\033[0m\n' "$*"; }

step "1/6 System packages"
sudo apt-get update
sudo apt-get install -y git build-essential cmake pkg-config wget usbutils \
    python3 python3-numpy python3-scipy \
    libusb-1.0-0-dev libncurses-dev libedit-dev libcurl4-openssl-dev \
    libfftw3-dev libmbedtls-dev libboost-program-options-dev libconfig++-dev libsctp-dev

step "2/6 GNU Radio and gr-gsm (GSM decoding)"
sudo apt-get install -y gnuradio gr-osmosdr || echo "!! GNU Radio packages failed - GSM tools will not work"
sudo apt-get install -y gr-gsm || echo "!! gr-gsm is not packaged for this Ubuntu - GSM tools need it built from source"

step "3/6 libbladeRF + bladeRF-cli (into /usr/local)"
if [ ! -x /usr/local/bin/bladeRF-cli ]; then
    rm -rf /tmp/bladeRF-src
    git clone --depth 1 --recursive https://github.com/Nuand/bladeRF.git /tmp/bladeRF-src
    mkdir -p /tmp/bladeRF-src/host/build && cd /tmp/bladeRF-src/host/build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DINSTALL_UDEV_RULES=ON \
          -DBLADERF_GROUP=plugdev ..
    make -j"$JOBS"
    sudo make install
    sudo ldconfig
    cd "$REPO"
else
    echo "already installed"
fi
sudo usermod -aG plugdev "$USER"
sudo udevadm control --reload-rules && sudo udevadm trigger || true

step "4/6 BladeRF 2.0 FPGA images (loaded automatically when the radio opens)"
FPGA_DIR=/usr/local/share/Nuand/bladeRF
sudo mkdir -p "$FPGA_DIR"
for img in hostedxA4 hostedxA9; do
    if [ ! -s "$FPGA_DIR/$img.rbf" ]; then
        sudo wget -q --show-progress -O "$FPGA_DIR/$img.rbf" "https://www.nuand.com/fpga/$img-latest.rbf" \
            || echo "!! could not download $img.rbf"
    fi
done

step "5/6 srsRAN 4G with our patch (LTE decoding)"
if [ ! -d "$HOME/srsRAN_4G" ]; then
    git clone --depth 1 https://github.com/srsran/srsRAN_4G.git "$HOME/srsRAN_4G"
fi
python3 "$REPO/osprey 2/srsran_watch/print_sib1.py"
mkdir -p "$HOME/srsRAN_4G/build" && cd "$HOME/srsRAN_4G/build"
cmake .. >/dev/null
make -j"$JOBS" pdsch_ue sib1_regen nr_decode nr_testgen
cd "$REPO"

step "6/6 Checks"
export LD_LIBRARY_PATH=/usr/local/lib
grep -q "LD_LIBRARY_PATH=/usr/local/lib" "$HOME/.bashrc" || echo 'export LD_LIBRARY_PATH=/usr/local/lib' >> "$HOME/.bashrc"
ls -l "$HOME/srsRAN_4G/build/lib/examples/pdsch_ue" && echo "srsRAN pdsch_ue: OK"
python3 lte_fastsync.py | tail -n 1
if lsusb | grep -qi "2cf0:5250"; then
    echo "BladeRF is plugged in:"; bladeRF-cli -p || true
else
    echo "BladeRF not plugged in (or not seen yet) - plug it into a USB 3 (blue) port"
fi
printf '\n\033[1;32mDone.\033[0m Log out and back in once (for the USB permission), then e.g.:\n'
echo "  cd $REPO && python3 lte_decode.py --earfcn 6300"
