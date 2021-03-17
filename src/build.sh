#!/bin/bash
if ! type git >/dev/null 2>&1; then
  sudo apt-get install git -y
fi
if [ ! -d "/home/ittc/wl8-build/" ]; then
  mkdir -p /home/ittc/wl8-build/ && cd /home/ittc/wl8-build/
  git clone git://git.ti.com/wilink8-wlan/build-utilites.git
fi
cd /home/ittc/wl8-build/build-utilites/
git checkout r8.8

export PATH=/opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/linux-devkit/sysroots/x86_64-arago-linux/usr/bin:$PATH

rm ./setup-env
sudo tee ./setup-env <<-'EOF'
export TOOLCHAIN_PATH=/opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/linux-devkit/sysroots/x86_64-arago-linux/usr/bin
export ROOTFS=./fs
export KERNEL_PATH=/opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/board-support/linux-4.19.94+gitAUTOINC+be5389fd85-gbe5389fd85
export CROSS_COMPILE=arm-linux-gnueabihf-
export ARCH=arm
[ "$TOOLCHAIN_PATH" != "" ] && export PATH=$TOOLCHAIN_PATH:$PATH
EOF
cp /opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/board-support/linux-4.19.94+gitAUTOINC+be5389fd85-gbe5389fd85/arch/arm/configs/tisdk_am335x-evm_defconfig /opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/board-support/linux-4.19.94+gitAUTOINC+be5389fd85-gbe5389fd85/arch/arm/configs/tisdk_am335x-evm_defconfigtest

sudo tee -a /opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/board-support/linux-4.19.94+gitAUTOINC+be5389fd85-gbe5389fd85/arch/arm/configs/tisdk_am335x-evm_defconfigtest <<-'EOF'
CONFIG_CFG80211=y
CONFIG_MAC80211=y
CONFIG_WLCORE=m
CONFIG_WLCORE_SDIO=m
CONFIG_WL18XX=m
CONFIG_NL80211_TESTMODE=y
CONFIG_MAC80211_MESH=y
CONFIG_MTD_RAM=y
CONFIG_SECURITY=y
CONFIG_MAC80211_DEBUGFS=y
CONFIG_CRYPTO_ARC4=y
CONFIG_CRYPTO_ECB=y
CONFIG_CRYPTO_MICHAEL_MIC=y
CONFIG_CRYPTO_CCM=y
CONFIG_CRYPTO_GCM=y
CONFIG_CRC7=y
CONFIG_INPUT_UINPUT=y

EOF

#./build_wl18xx.sh init

./verify_kernel_config.sh /opt/ti-processor-sdk-linux-am335x-evm-06.03.00.106/board-support/linux-4.19.94+gitAUTOINC+be5389fd85-gbe5389fd85/arch/arm/configs/tisdk_am335x-evm_defconfigtest
#./build_wl18xx.sh patch_kernel
#./build_wl18xx.sh kernel tisdk_am335x-evm_defconfig
#./build_wl18xx.sh update R8.8
