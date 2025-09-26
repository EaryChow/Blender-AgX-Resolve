# Who Wrote This Code?

This code was written by Jed Smith and Juan Pablo Zambrano. The base I believe was written by Jed Smith, and Juan Pablo Zambrano added the AgX mechanisms.

This repository includes some minor changes by Troy. And then the owner of this Blender-AgX fork included some changes for reproducing Blender's AgX result, also added more additional settings.

# AgX-Resolve

AgX Picture Formation for DaVinci Resolve

# Who

If you don't know what this is, it probably isn't for you just yet.

# What

Provides a flexible picture formation chain for DaVinci Resolve users.

# When

Currently testing out some things for folks who happen to find this repository. Not quite baked fully just yet.

# Where

In your DaVinci Resolve application `LUT` folder, create a subfolder called AgX. Place both files in there, and refresh your DaVinci Resolve installation.

# Why

You'll have to take it for a try.

# How

  1. Install as per Where above.
      1. Place the DCTL node. Select the Blender-AgX DCTL. Unlike the original AgX-Resolve repo, this fork defaults the input to Linear Rec.709, and then separately defines the working primaries, hence not requiring a pre-transform before the DCTL node.
2. Experiment.

Note that defaults are set for matching Blender's results.

That said, an author can set the values to whatever they choose, and the AgX mechanic will hold up.

# Parameters
These parameters work together, and thus one can expect any single parameter to drive other potential colour qualia.

## RGB Purity Attenuation
Greater values increase the rate of chromaticity attenuation as tristimulus values ascend. Lower values will slow the rate of attenuation. Slower rates may induce cognitively dissonant gradations, where a cognition mislocates the 2.5D "layer" from "under" to "over" or "through". For more information, please read the ramblings around [Picture Formation here](https://github.com/sobotka/scise/wiki/Picture-Formation).

## RGB Rotations
Controls the direction and rate of chromaticity angle flights. Value is degrees in CIE xy. Higher values will increase chromaticity angle flight speed toward the direction specified. 

## Inset
Purity Attenuation + RGB Rotation = Inset.

## Restore Purity
Controls the general purity of the primary when processing is complete. Some colourimetries will yield values that are too strong, breaking the surface of the picture. To reduce purity, decrease value. For a complete round trip no operation, set to match the Attenuation Rate above.

## Reverse RGB Rotation
Partially reverse the RGB Rotation, but due to the per-channel mechanism, it intentionally won't actually cancel the hue flight, but leave some "shifts" in place, which is the whole point of the Hue Flight mechanism.

## Outset
Restore Purity + Reverse RGB Rotation = Outset.

## Contrast
General contrasts of specified regions.

## Input Primaries
Input CIE colourimetry primaries for getting the values into the DCTL node

## Working Primaries
CIE colourimetry primaries for the working space. The input is transformed from the input primaries to the working primaries before all the AgX mechanism steps. This separates the input and the working space, making the node flexible for different input encodings.

## Input Transfer Characteristic Encoding
Input transfer characteristic that will be "undone" into linear before applying the Inset mechanism etc..

## Working Log Encoding
Log-like transfer for the working space.

## Output Primaries
Output primaries encoding.

## Contrast Pivot Offset
Offset added to the Working middle grey point in Log state. When offset is 0, the working Log middle grey is at where the working log would encode from linear 0.18, in othe words, lin_to_log(0.18).

## Log2 Min and Max Stop
These two settings are only used when the working Log is set to `Generic Log2`, which is a user controlled pure Log2 curve.

## Hue Flight Strength
The strength of the Per-Channel Hue Flight, 0 means completely chromaticity-linear, 1 means compeletely per-channel. Higher value would typically have yellower fire, for example. Blender's LUT sets it to 0.4.

## Tinting Scale and Tinting Hue
Mechanism ported from SB2383 python script. This adjusted the achromatic point * of the outset only * to create a tinting effect. It scales the achromatic point away from its original position, rotate it, then draw a line between the original position and the rotated position, and find the intersection point as "hull point". The final adjusted achromatic point derives from interpolating between the original position and the "hull point" position, with 0 meaning the original position, and 1 meaning the hull point. UI slider range is limited to prevent breaking of the image.

## Compensate for the Nagatives
This is the "Lower Guard Rail" used in Blender's AgX to handle negative values. 

## Log Encoded Output
Toggle to control whether to encode to the Working Log Encoding with Output Primaries for the output

## Use HDR
Checkbox to enable HDR-specific logic. When enabled, the middle grey of the image will be darkened. Darken by how much, depends on the HDR_Peak/SDR_Peak ratio determined by the following two settings.

## HDR Peak Nits
The peak nits for HDR output. UI limited to [400, 1000] with default of 1000.

## SDR Peak Nits
The SDR peak nits assumption. UI limted to [100, 400] with default of 203. Blender's LUT sets it to 100.

## HDR Extra Shoulder Power
Apply extra multiplier to sigmoid's shoulder power parameter, compensating for HDR's middle gray darkening.

## HDR purity
Adjust purity of the "HDR headroom". 

## HDR as percentage for SDR
In case `Use HDR` is checked with the output transfer function set to SDR ones. If checked, it outputs the image with a darkened middle grey. If not checked, it outputs the darkened-middle-grey image being scaled backup globally (which means clipping might occur). Default is checked to avoid output clipping. 
