import matplotlib.pyplot as plt
import numpy as np


class Braid:
    """
    Initialises the Braid class object, including calling the internal tracking
    function to generate strand positions.
    """

    def __init__(self, n, *ops):
        """
        Braid class object with internal strand tracking and drawing
        functionality.

        Attributes
        ----------
        braid_group : int
            The number of strands in the braid.
        braid_word : list
            Sequence of Artin operators that define the braid (Negative values
            indicate under-crossing strands).

        under-crossing_labels : list
            List of under-crossing strands' labels for each Artin operation
            respectively.
        """
        # Establish Group/No. of Strands
        self.braid_group = n

        # BraidWord and tracked strands developed at the same time
        self.braid_word = []
        for i in ops:
            # Check for invalid braid operation in selected group
            if abs(i) >= n:
                raise ValueError("Braid operations cannot exceed the containing braid group.")
            # Build braid word
            self.braid_word.append(i)

        self.undercrossing_labels = self._track_strands()

    def _track_strands(self):
        """
        Tracks strand positions through the given braid.

        top_labels : list
            List of 'braid_group' length, indicating strand names / labels
            present at the start, or 'top', of the braid. (Tracked internally
            from sequence of operations)
        bot_labels : list
            List of 'braid_group' length, indicating strand names / labels
            present at the end, or 'bottom', of the Braid. (Default
            [1, 2, ..., n])

        Returns
        -------
        under-crossing_labels : list
            List of under-crossing strands' labels for each Artin operation
            respectively.
        """
        # Initial Positions (AT BOTTOM)
        label_positions = [i + 1 for i in range(self.braid_group)]
        self.bot_labels = label_positions

        undercrossing_labels = []
        # Going through artin operators backwards (i.e. from bottom of braid)
        for op in reversed(self.braid_word):
            index = abs(op) - 1
            # Grab under-crossing string
            if np.sign(op) == -1:
                undercrossing_labels.append(label_positions[index + 1])
            else:
                undercrossing_labels.append(label_positions[index])
            # Swap
            (
                label_positions[index],
                label_positions[index + 1],
            ) = (
                label_positions[index + 1],
                label_positions[index],
            )

        # Final Positions (AT TOP)
        self.top_labels = label_positions

        undercrossing_labels.reverse()
        return undercrossing_labels

    def draw(self, style="comp", line_width=3, gap_size=3, color="rainbow", save=False):
        """
        Draws the Braid as a MatpLotlib figure.

        Parameters
        ----------
        style : "comp" or "ext"
            "comp" renders the image of the braid in a compact style with
            crossings parallel to one another if possible. "ext", for extended,
            shows the crossings in series.
        line_width : int (Default = 3)
            Thickness of the strands in the figure.
        gap_size : int (Default = 3)
            Amount of space shown at crossings for undercrossing strands.
        color : str
            Multicolor strands defined by "rainbow". Single fixed colour for
            all strands can be chosen from:
                {'b': blue,
                'g': green,
                'r': red,
                'c': cyan,
                'm': magenta,
                'y': yellow,
                'k': black,
                'w': white}
        """
        n = self.braid_group
        braid = self.braid_word

        fig, ax = plt.subplots(figsize=(4, 8))

        t = np.arange(0, 1.05, 0.05)

        # line break
        lb = np.array([[np.NaN, np.NaN]])

        x = 0
        layer = []

        strands = [np.empty((0, 2)) for _ in range(n)]

        # Start positions
        for s in range(1, n + 1):
            S = np.array([[s], [0]]).T
            strands[s - 1] = np.concatenate((strands[s - 1], S))

        # Calc crossing points
        for count, op in enumerate(braid):
            i = abs(op)

            if style == "comp":
                if i in layer or i + 1 in layer:
                    layer = []
                    x += 2
                layer.append(i)
                layer.append(i + 1)

            # Crossing Matrices
            R1 = np.array([0.5 * t**2 + i, -t - x]).T
            R2 = np.array([-0.5 * t**2 + 1 + i, t - 2 - x]).T
            L1 = np.array([-0.5 * t**2 + 1 + i, -t - x]).T
            L2 = np.array([0.5 * t**2 + i, t - 2 - x]).T

            if np.sign(op) == 1:  # over-crossing
                strands[i - 1] = np.concatenate((strands[i - 1], R1, np.flipud(R2)))
                strands[i] = np.concatenate(
                    (
                        strands[i],
                        L1[:-gap_size, :],
                        lb,
                        np.flipud(L2[:-gap_size, :]),
                    )
                )
            elif np.sign(op) == -1:  # under-crossing
                strands[i - 1] = np.concatenate(
                    (
                        strands[i - 1],
                        R1[:-gap_size, :],
                        lb,
                        np.flipud(R2[:-gap_size, :]),
                    )
                )
                strands[i] = np.concatenate((strands[i], L1, np.flipud(L2)))

            # swap lines for each operation
            strands[i - 1], strands[i] = (
                strands[i],
                strands[i - 1],
            )

            if style == "ext" or count + 1 == len(braid):
                x += 2

        # End positions
        for s in range(1, n + 1):
            S = np.array([[s], [-(x)]]).T
            strands[s - 1] = np.concatenate((strands[s - 1], S))
            if color == "rainbow":
                ax.plot(
                    strands[s - 1][:, 0],
                    strands[s - 1][:, 1],
                    lw=line_width,
                )
            else:
                ax.plot(
                    strands[s - 1][:, 0],
                    strands[s - 1][:, 1],
                    lw=line_width,
                    c=color,
                )

        # Figure Details
        fig.suptitle("Braid:   " + str(braid))
        ax.axis("off")
        ax.set_aspect(2 * n / x)
        plt.tight_layout()

        # SAVE CHECK
        if save:
            plt.savefig("test.svg")
        else:
            plt.show()

    def __str__(self):
        """
        Returns descriptive information about the Braid object.
        """
        return f"Braid: {self.braid_word}\nLabels: {self.undercrossing_labels}"
