#!/bin/bash

cd ../build
cmake .
make -j8
cd ../scripts
../build/bin/match_images_opencv_display  image_27.jpg image_26.jpg  match_27to26.png &
../build/bin/match_images_opencv_display  image_26.jpg image_27.jpg  match_26to27.png &
../build/bin/match_images_opencv_display  image_44.jpg image_45.jpg  match_44to45.png &
../build/bin/match_images_opencv_display  image_45.jpg image_44.jpg  match_45to44.png &
